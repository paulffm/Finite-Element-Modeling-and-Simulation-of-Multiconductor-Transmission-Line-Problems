# coding=utf-8
"""File containing the implementation of DistributedSolidConductor

.. sectionauthor:: Bundschuh
"""

from typing import Tuple, TYPE_CHECKING, Union, List, Dict
from copy import deepcopy

import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix

from pyrit import get_logger

from pyrit.region import Regions
from . import FieldCircuitCoupling, SolidConductor, Excitations

if TYPE_CHECKING:
    from pyrit.excitation.FieldCircuitCoupling import ShrinkInflateProblem

ComplexNumber = Union[int, float, complex]
Number = Union[int, float]

logger = get_logger(__name__)


class DistributedSolidConductor(FieldCircuitCoupling):
    """A set of solid conductors that belong to the same component (from an outside view)

    An example is a wire coil where each wire is resolved in the geometry and by the mesh. Then each wire is modelled as
    a solid conductor. However, all the single wires belong to one component (a wire coil in this case).

    *How to use this conductor model?*

    It is mandatory to define a number of single physical groups (one for each conductor inside
    DistributedSolidConductor). Assign all these physical groups to one instance of DistributedSolidConductor and make
    sure that the number of physical groups coincides with the `number_conductors` from the constructor.

    >>> turns = 50
    >>> geo = Geometry('model')
    >>> conductors = [geo.create_physical_group(100 + k, 2, f'conductor_{k}') for k in range(turns)]
    >>> dsc = DistributedSolidConductor(turns, current=1+1j, name='coil')
    >>> geo.add_excitation_to_physical_group(dsc, *conductors)

    If you use the pyrit to build the geometry, you have to add an entity to each physical group.
    """

    def __init__(self, number_conductors: int, *, current=None, voltage=None, name=''):
        super().__init__(current, voltage, name)

        self.number_conductors = number_conductors
        self.solid_conductors: List[SolidConductor] = []

    @property
    def additional_equations(self) -> int:
        return 2 * self.number_conductors

    def _get_region_id(self, regions):
        region_ids = super()._get_region_id(regions)

        if len(region_ids) != self.number_conductors:
            logger.warning("The number of regions (%d) does not match the number of conductors (%d). The number of"
                           "conductors is updated.", len(region_ids), self.number_conductors)
            self.number_conductors = len(region_ids)

        return region_ids

    def set_solution(self, solution: np.ndarray):
        """Get the solution for the distributed solid conductors and processes it

        Parameters
        ----------
        solution : np.ndarray
            An array of length 2 times the number of conductors where the voltage is in the odd and the current in the
            even elements.
        """
        if not len(solution) == self.additional_equations:
            raise ValueError("The length of solution is not equal to the number of additional equations")

        for k in range(0, self.number_conductors):
            self.solid_conductors[k].set_solution(solution[2 * k:2 * k + 2])

        if self.voltage_mode:
            self._current = sum((sc.current for sc in self.solid_conductors))
        if self.current_mode:
            self._voltage = sum((sc.voltage for sc in self.solid_conductors))

    def single_currents(self) -> Dict[int, complex]:
        """Return a dict with one current per solid conductor.

        Returns
        -------
        single_currents : Dict[int, complex]
            Key is the id of the region and value is the current through the solid conductor defined on this region
        """
        # noinspection PyUnresolvedReferences
        return {sc.region_id: sc.current for sc in self.solid_conductors}

    def single_voltages(self) -> Dict[int, complex]:
        """Return a dict with one voltage per solid conductor.

        Returns
        -------
        single_voltages : Dict[int, complex]
            Key is the id of the region and value is the voltage over the solid conductor defined on this region
        """
        # noinspection PyUnresolvedReferences
        return {sc.region_id: sc.voltage for sc in self.solid_conductors}

    def single_impedances(self) -> Dict[int, complex]:
        """Return a dict with one impedance per solid conductor.

        Returns
        -------
        single_impedances : Dict[int, complex]
            Key is the id of the region and value is the impedance of the solid conductor defined on this region
        """
        # noinspection PyUnresolvedReferences
        return {sc.region_id: sc.impedance for sc in self.solid_conductors}

    def single_resistances(self) -> Dict[int, float]:
        """Return a dict with one resistance per solid conductor.

        Returns
        -------
        single_resistances : Dict[int, complex]
            Key is the id of the region and value is the resistance of the solid conductor defined on this region
        """
        # noinspection PyUnresolvedReferences
        return {sc.region_id: sc.resistance for sc in self.solid_conductors}

    def single_inductances(self, omega: float) -> Dict[int, float]:
        """Return a dict with one inductance per solid conductor.

        Returns
        -------
        single_inductances : Dict[int, complex]
            Key is the id of the region and value is the inductance of the solid conductor defined on this region
        """
        # noinspection PyUnresolvedReferences
        return {sc.region_id: sc.inductance(omega) for sc in self.solid_conductors}

    def get_coupling_values(self, problem: 'ShrinkInflateProblem') -> \
            Tuple[coo_matrix, coo_matrix, coo_matrix, coo_matrix]:
        region_ids = self._get_region_id(problem.regions)

        tmp_problem = deepcopy(problem)

        bottom_matrices, right_matrices, diagonal_matrices, rhs_matrices = [], [], [], []
        self.solid_conductors = []
        for k, region_id in enumerate(region_ids):
            sc = SolidConductor(current=self.current, voltage=self.voltage,
                                name='_'.join([self.name, f"conductor_{k}"]))
            sc.region_id = region_id
            self.solid_conductors.append(sc)

            region = deepcopy(problem.regions.get_regi(region_id))
            region.exci = sc.ID

            excitations_tmp = Excitations(sc)
            regions_tmp = Regions(region)

            tmp_problem.regions = regions_tmp
            tmp_problem.excitations = excitations_tmp

            matrix_bot_tmp, matrix_right_tmp, matrix_diag_tmp, matrix_rhs_tmp = sc.get_coupling_values(tmp_problem)
            sc.cleanup(tmp_problem)

            bottom_matrices.append(matrix_bot_tmp)
            right_matrices.append(matrix_right_tmp)
            diagonal_matrices.append(matrix_diag_tmp)
            rhs_matrices.append(matrix_rhs_tmp)

        bottom_matrix: sparse.csr_matrix = sparse.vstack(bottom_matrices, format='coo')
        right_matrix: sparse.csr_matrix = sparse.hstack(right_matrices, format='coo')
        diagonal_matrix = sparse.block_diag(diagonal_matrices, format='coo')
        rhs_fcc = sparse.vstack(rhs_matrices, format='coo')

        return bottom_matrix, right_matrix, diagonal_matrix, rhs_fcc

    def cleanup(self, problem: 'ShrinkInflateProblem') -> None:
        pass
