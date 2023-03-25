# coding=utf-8
"""File containing the class FieldCircuitCoupling

.. sectionauthor:: Bundschuh
"""

from abc import abstractmethod
from typing import TYPE_CHECKING, Protocol, Tuple, Dict, Generator, Union, Any
from dataclasses import dataclass
import numpy as np
from numpy import ndarray
from scipy.sparse import coo_matrix
from scipy import sparse

from pyrit.excitation.Exci import Exci

from pyrit import get_logger

if TYPE_CHECKING:
    from pyrit.mesh import Mesh
    from pyrit.shapefunction import EdgeShapeFunction
    from pyrit.region import Regions
    from pyrit.bdrycond import BdryCond
    from pyrit.material import Materials

    from . import Excitations

ComplexNumber = Union[int, float, complex]
Number = Union[int, float]

logger = get_logger(__name__)


@dataclass
class ShrinkInflateProblem(Protocol):
    """Protocol for a problem that meets the requirements for shrinking and inflating."""

    mesh: 'Mesh'
    shape_function: 'EdgeShapeFunction'
    regions: 'Regions'
    materials: 'Materials'
    boundary_conditions: 'BdryCond'
    excitations: 'Excitations'
    frequency: float
    omega: float


class FieldCircuitCoupling(Exci):
    r"""Abstract class for field-circuit coupling.

    Field-circuit coupling means to couple a FE model to an external electrical circuit. This class does this coupling.
    At the time it only works to have a single voltage or current source attached to the FE modelled device. For that,
    you can specify the current or the voltage, that is delivered by the source, in the constructor.

    It is important to use shrink and inflate of this class in order to build the matrices appropriately. When using
    the shrink and inflate method of the boundary conditions, this will lead to errors.

    Examples
    --------
    We suppose that we have an instance of an appropriate problem class ``problem``. In the excitation object stored in
    ``problem`` there is a field-circuit coupling excitation. Suppose it was created from a class ``FCC`` that inherits
    from ``FieldCircuitCoupling``. Its creation could look like this:

    >>> fcc = FCC(current=5)

    So the current through the device is set to :math:`5\,\mathrm{A}`. With ``matrix`` and ``rhs`` being the matrix and
    the right-hand side vecotr from a shape function object, the problem is shrunk and inflated by

    >>> matrix_shrink, rhs_shrink = FieldCircuitCoupling.shrink(matrix,rhs,problem)
    >>> ...
    >>> solution_vector, fcc_ids = FieldCircuitCoupling.inflate(solution_shrink, problem)

    The ``solution_vector`` then is the solution of the original system of equations and ``fcc_ids`` is a list of
    indices of field-circuit elements. In this case the list would contain only th ID of ``fcc``. The solution if
    written to this object (in this case ``fcc``). So if the current was given, the voltage is written to the object
    and vice versa.
    """

    def __init__(self, current: ComplexNumber = None, voltage: ComplexNumber = None,
                 name: str = ''):
        """Initializer of FieldCircuitCoupling.

        Either a current or  a voltage has to be given.

        Parameters
        ----------
        current : ComplexNumber, optional
            The current through the component. Default is None
        voltage : ComplexNumber, optional
            The voltage over the component. Default is None
        name : str, optional
            The name of the component. Default is ''
        """
        super().__init__(None, name=name)
        self._mode = None
        self._check_current_voltage(current, voltage)
        self._current = current
        self._voltage = voltage

    def _check_current_voltage(self, current: Any, voltage: Any):
        """Check if either a current or a voltage is given. Set the mode appropriately.

        Parameters
        ----------
        current : Any
            A current or None
        voltage : Any
            A voltage or None
        """
        if (current is None and voltage is None) or (current is not None and voltage is not None):
            raise ValueError("You have to give either a current or a voltage")
        self._mode = 'current' if voltage is None else 'voltage'

    @property
    @abstractmethod
    def additional_equations(self) -> int:
        """Returns the number of equations added to the standard finite element model."""

    def _get_region_id(self, regions):
        """Get the region of the own excitation."""
        return regions.find_regions_of_excitations(self.ID)

    @abstractmethod
    def set_solution(self, solution: np.ndarray):
        """Get the solution for the component and processes it

        Parameters
        ----------
        solution : np.ndarray
            An array of length 2.
        """

    def activate_voltage_mode(self):
        """Switch to voltage mode"""
        self._mode = 'voltage'

    def activate_current_mode(self):
        """Switch to current mode"""
        self._mode = 'current'

    @property
    def voltage_mode(self):
        """True if the component is in voltage mode, i.e. the voltage is given."""
        return self._mode == 'voltage'

    @property
    def current_mode(self):
        """True if the component is in current mode, i.e. the current is given."""
        return self._mode == 'current'

    @property
    def current(self):
        """The current through the component."""
        return self._current

    @current.setter
    def current(self, current):
        if not self.current_mode:
            logger.warning('The current was set while in voltage mode.')
        self._current = current

    @property
    def voltage(self):
        """The voltage over the component."""
        return self._voltage

    @voltage.setter
    def voltage(self, voltage):
        if not self.voltage_mode:
            logger.warning('The voltage was set while in current mode.')
        self._voltage = voltage

    @property
    def impedance(self):
        """The impedance of the component."""
        return self.voltage / self.current

    @property
    def resistance(self):
        """The resistance of the component."""
        return self.impedance.real

    def inductance(self, omega: float):
        """The inductance of the component."""
        return self.impedance.imag / omega

    @abstractmethod
    def get_coupling_values(self, problem: ShrinkInflateProblem) -> \
            Tuple[coo_matrix, coo_matrix, coo_matrix, coo_matrix]:
        r"""Return all matrices and values needed for the field circuit coupling.

        Let :math:`\mathbf{A}\mathbf{x} = \mathbf{b}` the system to shrink. Let furthermore be :math:`\mathbf{B}` the
        *bottom matrix*, :math:`\mathbf{R}` the *right matrix*, :math:`c` the *coupling value* and :math:`r` the
        *right-hand-side value*. Then the new matrix is of the form

        .. math::

            \begin{pmatrix} \mathbf{A} & \mathbf{R}\\\mathbf{B} & c\end{pmatrix}
            \begin{pmatrix}\mathbf{x}\\\xi\end{pmatrix}=\begin{pmatrix}\mathbf{b}\\r\end{pmatrix}

        Parameters
        ----------
        problem : ShrinkInflateProblem

        Returns
        -------
        bottom_matrix : coo_matrix
            The matrix that goes below the original matrix.
        right_matrix : coo_matrix
            The matrix that goes right of the original matrix.
        coupling value : coo_matrix
            The matrix that goes in the bottom right corner.
        rhs_value : coo_matrix
            The matrix for the right-hand-side.
        """

    # region Exterior functionality

    @classmethod
    def iter_fcc(cls, excitations: 'Excitations') -> Generator['FieldCircuitCoupling', None, None]:
        """A generator for iterating over all Excis in Excitations of the same class.

        Parameters
        ----------
        excitations : Excitations
            An excitations object.

        Yields
        ------
        fcc : cls
            An instance of the current class.
        """
        for fcc in excitations:
            if isinstance(fcc, cls):
                yield fcc

    @classmethod
    def get_fccs(cls, excitations: 'Excitations') -> Dict[int, 'FieldCircuitCoupling']:
        """Get a dictionary with instances of the class

        Parameters
        ----------
        excitations : 'Excitations'
            An excitations object

        Returns
        -------
        fccs : Dict[int,'FieldCircuitCoupling']
        """
        return {fcc.ID: fcc for fcc in cls.iter_fcc(excitations)}

    @classmethod
    def shrink(cls, matrix: coo_matrix, rhs: coo_matrix, problem: ShrinkInflateProblem, integration_order: int = 1) -> \
            Tuple[coo_matrix, coo_matrix]:
        """Shrink a system of equations with respect to boundary conditions and field circuit coupling.

        First the system in shrunk as normal (as you would do it with *BdryCond*). The resulting system of equations is
        then further modified with respect to the field circuit coupling.

        As a consequence, if you want to perform a field circuit coupling, you only need to use this shrink method
        instead of the one from *BdryCond*.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the original system of equations.
        rhs : coo_matrix
            The right-hand-side of the original system of equations.
        problem : ShrinkInflateProblem
            A problem object
        integration_order : int, optional
            The integration order used in subsequent methods.

        Returns
        -------
        matrix_shrink : coo_matrix
            The shrunk matrix.
        rhs_shrink : coo_matrix
            The shrunk right-hand-side.
        """
        # compute the coupling values for every fcc
        bottom_matrices, right_matrices, diagonal_matrices, rhs_matrices = [], [], [], []

        try:
            for fcc in cls.iter_fcc(problem.excitations):
                # get the coupling values
                # noinspection PyTupleAssignmentBalance
                matrix_bot_tmp, matrix_right_tmp, matrix_diag_tmp, matrix_rhs_tmp = fcc.get_coupling_values(problem)
                bottom_matrices.append(matrix_bot_tmp)
                right_matrices.append(matrix_right_tmp)
                diagonal_matrices.append(matrix_diag_tmp)
                rhs_matrices.append(matrix_rhs_tmp)
        finally:  # Clean up
            for fcc in cls.iter_fcc(problem.excitations):
                fcc.cleanup(problem)

        # Assembly of the matrices
        bottom_matrix: sparse.csr_matrix = sparse.vstack(bottom_matrices, format='csr')
        right_matrix: sparse.csr_matrix = sparse.hstack(right_matrices, format='csr')
        diagonal_matrix = sparse.block_diag(diagonal_matrices, format='csr')
        rhs_fcc = sparse.vstack(rhs_matrices, format='csr')

        # shrink the system with the sft
        ret = problem.shape_function.shrink(matrix, rhs, problem, integration_order=integration_order)
        matrix_shrink_sf, rhs_shrink_sf, indices_not_dof, _, _ = ret

        # compose the coupling matrices and coupling values with the information from the sft.shrink
        all_indices = np.arange(rhs.shape[0])
        # num_fcc = len(bottom_matrices)
        new_lines = bottom_matrix.shape[0]
        new_shape = matrix_shrink_sf.shape[0]
        indices_dof = np.setdiff1d(all_indices, indices_not_dof)
        bottom_matrix = bottom_matrix[np.ix_(np.arange(new_lines), indices_dof)]
        right_matrix = right_matrix[np.ix_(indices_dof, np.arange(new_lines))]
        bottom_matrix.resize((bottom_matrix.shape[0], new_shape))
        right_matrix.resize((new_shape, right_matrix.shape[1]))

        # compute the final system
        matrix_shrink = sparse.bmat([[matrix_shrink_sf, right_matrix],
                                     [bottom_matrix, diagonal_matrix]], format='coo')
        rhs_shrink = sparse.vstack([rhs_shrink_sf, rhs_fcc], format='coo')

        return matrix_shrink, rhs_shrink

    @classmethod
    def inflate(cls, solution: ndarray, problem: ShrinkInflateProblem) -> Tuple[ndarray, Dict[int, ComplexNumber]]:
        """Inflate the solution.

        After the shrunk system is solved, use this function for inflating the solution. This contains inflation with
        respect to boundary conditions and field circuit coupling.

        Parameters
        ----------
        solution : ndarray
            The solution of the shrunk system.
        problem : ShrinkInflateProblem
            A problem object.

        Returns
        -------
        solution_inflated : ndarray
            The inflated solution without any artifacts from field circuit coupling.
        solution_values : List[int]
            A List with the IDs of the FieldCircuitCoupling instances. The solutions are written to those instances.
        """
        # Get a List of IDs of the fcc
        ids = []  # List of all IDs of the fccs
        additional_equations = []  # List of the number of additional equations per fcc

        # Determine the two lists
        for fcc in cls.iter_fcc(problem.excitations):
            ids.append(fcc.ID)
            additional_equations.append(fcc.additional_equations)

        num_additional_equations = sum(additional_equations)
        solutions_fcc = solution[-num_additional_equations:]  # Vector of the fcc solutions

        # Assign the fcc solution to the fcc
        start = 0
        for k, fcc in enumerate(cls.iter_fcc(problem.excitations)):
            solution_fcc = solutions_fcc[start:start + additional_equations[k]]
            start += additional_equations[k]
            fcc.set_solution(solution_fcc)

        # Standard inflate
        solution_inflated = problem.shape_function.inflate(solution[:-num_additional_equations], problem)

        return solution_inflated, ids

    @abstractmethod
    def cleanup(self, problem: ShrinkInflateProblem) -> None:
        """Delete the added temporary fields in problem.

        Parameters
        ----------
        problem : ShrinkInflateProblem
            A problem object.
        """

    # endregion
