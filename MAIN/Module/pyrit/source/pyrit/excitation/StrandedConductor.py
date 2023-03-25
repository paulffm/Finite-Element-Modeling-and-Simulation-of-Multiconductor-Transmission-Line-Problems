# coding=utf-8
"""Implementation of StrandedConductor.

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING, Tuple, Union

import numpy as np
from scipy.sparse import coo_matrix

from pyrit import get_logger

from pyrit.material import Resistivity
from . import FieldCircuitCoupling

logger = get_logger(__name__)
if TYPE_CHECKING:
    from .FieldCircuitCoupling import ShrinkInflateProblem

ComplexNumber = Union[int, float, complex]
Number = Union[int, float]


class StrandedConductor(FieldCircuitCoupling):
    r"""Class representing a stranded conductor.

    In a stranded conductor, the current density is constant even for elevated frequencies. This is due to the
    assumption that the conductor consists of many single strands that are insulated to each other.

    Notes
    -----
    The current distribution function :math:`\vec{\chi}` used in this class has to satisfy

    .. math::

        \int_{S_c} \vec{\chi}\cdot\mathrm{d}\vec{S} = 1\,,

    with :math:`S_c` being the cross-section of the conducting domain. Note the 1 on the right-hand side. This may be
    defined differently. One could also define :math:`\vec{\chi}` such that there stands the number of windings
    :math:`N` on the right-hand side. Since that is here not the case, the current density is

    .. math::

        \vec{J} = NI\vec{\chi}\,.
    """

    # Names of the attributes that are temporarily saved in problem
    __mass_one_name = 'fcc_stranded_conductor_mass_one'
    __mass_resistivity_name = 'fcc_stranded_conductor_mass_resistivity'

    def __init__(self, *, current: ComplexNumber = None, voltage: ComplexNumber = None, windings: int = 10,
                 name: str = None):
        """Constructor of a stranded conductor.

        Either a current or  a voltage has to be given.

        Parameters
        ----------
        current : ComplexNumber, optional
            The current through the component. Default is None
        voltage : ComplexNumber, optional
            The voltage over the component. Default is None
        windings : int, optional
            The number of windings/strands in the conductor. Default is 1
        name : str, optional
            The name of the conductor.
        """
        super().__init__(current, voltage, name)
        self.windings = windings
        self.current_distribution_function = None

        self.base_coupling_matrix: coo_matrix = None
        self.bottom_coupling_matrix: coo_matrix = None
        self.right_coupling_matrix: coo_matrix = None
        self.coupling_value = None

    @property
    def additional_equations(self) -> int:
        return 2

    def set_solution(self, solution: np.ndarray):
        """Get the solution for the stranded conductor and processes it

        Parameters
        ----------
        solution : np.ndarray
            An array of length 2 where the current is in the first element and the voltage in the second.
        """
        if not len(solution) == 2:
            raise ValueError("There must be 2 solutions.")
        current, voltage = solution[0], solution[1]
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
            diagonal_matrix = coo_matrix(np.array([[self.windings ** 2 * self.coupling_value, -1], [0, 1]]))
            rhs_matrix = coo_matrix(np.array([[0], [self.voltage]]))
        else:
            diagonal_matrix = coo_matrix(np.array([[self.windings ** 2 * self.coupling_value, -1], [1, 0]]))
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

        # Get the one-mass matrix. Get it from problem or compute it.
        if hasattr(problem, self.__mass_one_name):
            mass_one = problem.__getattribute__(self.__mass_one_name)
        else:
            mass_one = problem.shape_function.mass_matrix(1)
            problem.__setattr__(self.__mass_one_name, mass_one)
        # Get the mass matrix. Get it form problem or compute it.
        if hasattr(problem, self.__mass_resistivity_name):
            mass = problem.__getattribute__(self.__mass_resistivity_name)
        else:
            mass = problem.shape_function.mass_matrix(problem.regions, problem.materials, Resistivity)
            problem.__setattr__(self.__mass_resistivity_name, mass)

        # Compute the current distribution function if it has not been computed yet.
        if self.current_distribution_function is None:
            self.compute_current_distribution_function(problem)

        self.compute_coupling_matrix(mass_one, problem.omega)
        self.compute_coupling_value(mass)

    def compute_coupling_matrix(self, mass_one: coo_matrix, omega: float) -> None:
        """Compute the coupling matrix and save it internally.

        Parameters
        ----------
        mass_one : coo_matrix
            The one-mass matrix.
        omega : float
            The angular frequency.
        """
        if self.base_coupling_matrix is None:
            self.base_coupling_matrix = mass_one @ self.current_distribution_function
            self.base_coupling_matrix = self.base_coupling_matrix.tocoo()
        self.bottom_coupling_matrix = 1j * omega * self.windings * self.base_coupling_matrix.transpose()
        self.right_coupling_matrix = -1 * self.windings * self.base_coupling_matrix

    def compute_coupling_value(self, mass):
        """Compute the coupling value und save it internally.

        Parameters
        ----------
        mass : coo_matrix
            The mass matrix.
        """
        coupling_value = self.current_distribution_function.transpose() @ mass @ self.current_distribution_function
        self.coupling_value = coupling_value.data[0]

    def compute_current_distribution_function(self, problem: 'ShrinkInflateProblem'):
        """Compute the current distribution function and save it internally.

        Parameters
        ----------
        problem : pyrit.excitation.FieldCircuitCoupling.ShrinkInflateProblem
            A problem object.
        """
        self.current_distribution_function = problem.shape_function.current_distribution_function(problem.regions, self)

    def cleanup(self, problem: 'ShrinkInflateProblem'):
        self.base_coupling_matrix = None
        self.bottom_coupling_matrix = None
        self.right_coupling_matrix = None
        self.coupling_value = None
        self.current_distribution_function = None

        if hasattr(problem, self.__mass_one_name):
            problem.__delattr__(self.__mass_one_name)
        if hasattr(problem, self.__mass_resistivity_name):
            problem.__delattr__(self.__mass_resistivity_name)
