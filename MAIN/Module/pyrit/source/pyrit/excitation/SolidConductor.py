# coding=utf-8
"""Implementation of SolidConductor.

.. sectionauthor:: Bundschuh
"""

from typing import Tuple, TYPE_CHECKING, Union

import numpy as np
from scipy.sparse import coo_matrix

from pyrit import get_logger

from pyrit.material import Conductivity
from . import FieldCircuitCoupling

if TYPE_CHECKING:
    from pyrit.excitation.FieldCircuitCoupling import ShrinkInflateProblem

ComplexNumber = Union[int, float, complex]
Number = Union[int, float]

logger = get_logger(__name__)


class SolidConductor(FieldCircuitCoupling):
    """Class representing a solid conductor.

    The current density is allowed to redistribute inside a solid conductor.
    """

    # Name of the attribute that is temporarily saved in problem
    __mass_name = 'fcc_solid_conductor_mass_conductivity'

    def __init__(self, *, current: ComplexNumber = None, voltage: ComplexNumber = None, name: str = None):
        """Constructor of a solid conductor.

        Either a current or  a voltage has to be given.

        Parameters
        ----------
        current : ComplexNumber, optional
            The current through the component. Default is None
        voltage : ComplexNumber, optional
            The voltage over the component. Default is None
        name : str, optinal
            The name of the solid conductor. Default is None.
        """
        super().__init__(current, voltage, name)
        self.voltage_distribution_function = None

        self.base_coupling_matrix: coo_matrix = None
        self.bottom_coupling_matrix: coo_matrix = None
        self.right_coupling_matrix: coo_matrix = None
        self.coupling_value = None

    @property
    def additional_equations(self) -> int:
        return 2

    def set_solution(self, solution: np.ndarray):
        """Get the solution for the solid conductor and processes it

        Parameters
        ----------
        solution : np.ndarray
            An array of length 2 where the voltage is in the first element and the current in the second.
        """
        if not len(solution) == 2:
            raise ValueError("There must be 2 solutions.")
        voltage, current = solution[0], solution[1]
        if self.voltage_mode:
            if not np.isclose(voltage, self.voltage):
                logger.warning("The voltage was given but has changed in the solution vector")
            self._current = current

        if self.current_mode:
            if not np.isclose(current, self.current):
                logger.warning("The current was given but has changed in the solution vector")
            self._voltage = voltage

    def get_coupling_values(self, problem: 'ShrinkInflateProblem') -> \
            Tuple[coo_matrix, coo_matrix, coo_matrix, coo_matrix]:
        if self.bottom_coupling_matrix is None or self.right_coupling_matrix is None or self.coupling_value is None:
            self.compute_coupling_values(problem)

        bottom_matrix: coo_matrix = self.bottom_coupling_matrix.copy()
        bottom_matrix.resize((2, self.bottom_coupling_matrix.shape[1]))

        right_matrix: coo_matrix = self.right_coupling_matrix.copy()
        right_matrix.resize((self.right_coupling_matrix.shape[0], 2))

        if self.voltage_mode:
            diagonal_matrix = coo_matrix(np.array([[self.coupling_value, -1], [1, 0]]))
            rhs_matrix = coo_matrix(np.array([[0], [self.voltage]]))
        else:
            diagonal_matrix = coo_matrix(np.array([[self.coupling_value, -1], [0, 1]]))
            rhs_matrix = coo_matrix(np.array([[0], [self.current]]))

        return bottom_matrix, right_matrix, diagonal_matrix, rhs_matrix

    def compute_coupling_values(self, problem: 'ShrinkInflateProblem'):
        """Computes matrices needed for the coupling and stored them internally.

        Parameters
        ----------
        problem : ShrinkInflateProblem
            A problem instance.
        """
        if len(self._get_region_id(problem.regions)) != 1:
            raise Exception("This conductor is defined on more than one region. This does not make sense.")

        if hasattr(problem, self.__mass_name):
            mass = problem.__getattribute__(self.__mass_name)
        else:
            mass = problem.shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)
            problem.__setattr__(self.__mass_name, mass)

        if self.voltage_distribution_function is None:
            self.compute_voltage_distribution_function(problem)
        self.compute_coupling_matrix(mass, problem.omega)
        self.compute_coupling_value(mass)

    def compute_coupling_matrix(self, mass: coo_matrix, omega: float) -> None:
        """Compute the coupling matrix and save it internally.

        Parameters
        ----------
        mass : coo_matrix
            The mass matrix.
        omega : float
            Angular frequency.
        """
        if self.base_coupling_matrix is None:
            self.base_coupling_matrix = mass @ self.voltage_distribution_function
            self.base_coupling_matrix = self.base_coupling_matrix.tocoo()
        self.bottom_coupling_matrix = -1j * omega * self.base_coupling_matrix.transpose()
        self.right_coupling_matrix = -1 * self.base_coupling_matrix

    def compute_coupling_value(self, mass: coo_matrix) -> None:
        """Compute the coupling value und save it internally.

        Parameters
        ----------
        mass : coo_matrix
            The mass matrix.
        """
        if self.right_coupling_matrix is None:
            coupling_value = self.voltage_distribution_function.transpose() @ mass @ self.voltage_distribution_function
        else:
            coupling_value = -1 * self.voltage_distribution_function.transpose() @ self.right_coupling_matrix
        self.coupling_value = coupling_value.data[0]

    def compute_voltage_distribution_function(self, problem: 'ShrinkInflateProblem') -> None:
        """Compute the voltage distribution function and save it internally.

        Parameters
        ----------
        problem : pyrit.excitation.FieldCircuitCoupling.ShrinkInflateProblem
            A problem object.
        """
        self.voltage_distribution_function = problem.shape_function.voltage_distribution_function(problem.regions, self)

    def cleanup(self, problem: 'ShrinkInflateProblem'):
        self.base_coupling_matrix = None
        self.bottom_coupling_matrix = None
        self.right_coupling_matrix = None
        self.coupling_value = None
        self.voltage_distribution_function = None

        if hasattr(problem, self.__mass_name):
            problem.__delattr__(self.__mass_name)
