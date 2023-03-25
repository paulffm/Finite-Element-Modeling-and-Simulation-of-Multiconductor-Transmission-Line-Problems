# coding=utf-8
"""
File containing the boundary condition class.

<<<<<<< HEAD
.. sectionauthor:: Bundschuh
=======
.. sectionauthor:: bundschuh
>>>>>>> PYR-94
"""

from typing import List, Tuple, Dict, TYPE_CHECKING, Any
from dataclasses import dataclass

from numpy import ndarray
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix

from . import BCBinary, BCDirichlet, BCFloating, BCNeumann, BCRobin

if TYPE_CHECKING:
    from pyrit.bdrycond import BC
    from pyrit.mesh import Mesh
    from pyrit.region import Regions
    from pyrit.shapefunction import ShapeFunction


class BdryCond:
    """Container class for collecting different boundary conditions for a problem."""

    def __init__(self, *bcs: 'BC'):
        """
        Constructor of BdryCond

        Parameters
        ----------
        bcs : BC, optional
            A number of instances of BC
        """
        self.__bcs = {}
        if len(bcs) > 0:
            self.add_bc(*bcs)

    def __iter__(self):
        return self.__bcs.values().__iter__()

    def add_bc(self, *bcs: 'BC') -> None:
        """
        Adding BCs to BdryCond.

        Parameters
        ----------
        bcs : BC
            A number of instances of BC

        Returns
        -------
        None
        """
        for bc in bcs:
            self.__bcs[bc.ID] = bc

    def delete_bc(self, *IDs: int) -> None:
        """
        Delete a number of BCs from BdryCond

        Parameters
        ----------
        IDs : int
            IDs of the BCs to delete

        Returns
        -------
        None
        """
        for ID in IDs:
            del self.__bcs[ID]

    def get_bc(self, ID: int) -> 'BC':
        """
        Return the BC specified by ID

        Parameters
        ----------
        ID : int
            The ID of the BC

        Returns
        -------
        bc: BC
            A instance of BC
        """
        return self.__bcs[ID]

    def get_ids(self) -> List[int]:
        """
        Returns a list of all IDs

        Returns
        -------
        IDs: List[int]
            The List
        """
        return list(self.__bcs.keys())

    def num_bc(self) -> int:
        """
        Returns the number of stored BCs.

        Returns
        -------
        num_bc: int
            The number of stored BCs
        """
        return len(self.get_ids())

    @property
    def is_linear(self) -> bool:
        """Returns True if all boundary conditions are linear.

        Returns
        -------
        linear : bool
        """
        return all((bc.is_linear for bc in self.__bcs.values()))

    @property
    def is_time_dependent(self) -> bool:
        """Returns True if at least one boundary condition is time dependent.

        Returns
        -------
        time_dependent : bool
        """
        return any((bc.is_time_dependent for bc in self.__bcs.values()))

    @property
    def is_homogeneous(self) -> bool:
        """Returns True if all boundary conditions are homogeneous.

        Returns
        -------
        homogeneous : bool
        """
        return all((bc.is_homogeneous for bc in self.__bcs.values()))

    @property
    def is_constant(self) -> bool:
        """Returns True if all boundary conditions are constant.

        Returns
        -------
        constant : bool
        """
        return all((bc.is_constant for bc in self.__bcs.values()))

    def update(self, name: str, val: Any):
        """Update a value to a name.

        Calls the function :py:meth:`pyrit.bdrycond.BC.update` on every boundary condition.

        Parameters
        ----------
        name : str
            The name of the field to update.
        val : Any
            The value to update.
        """
        for bc in self.__bcs.values():
            bc.update(name, val)

    def dict_of_boundary_condition(self) -> Dict[str, List[int]]:
        """Returns a dict with the IDs of the boundary conditions, divided by type.

        Returns
        -------
        dict_of_bc : Dict[str, List[int]]
            Dict with boundary conditions.
        """
        list_dirichlet, list_neumann, list_robin, list_floating, list_binary = [], [], [], [], []
        for bc in self.__bcs.values():
            if isinstance(bc, BCDirichlet):
                list_dirichlet.append(bc.ID)
            elif isinstance(bc, BCNeumann):
                list_neumann.append(bc.ID)
            elif isinstance(bc, BCRobin):
                list_robin.append(bc.ID)
            elif isinstance(bc, BCFloating):
                list_floating.append(bc.ID)
            elif isinstance(bc, BCBinary):
                list_binary.append(bc.ID)
        return {'dirichlet': list_dirichlet, 'neumann': list_neumann, 'robin': list_robin, 'floating': list_floating,
                'binary': list_binary}

    def regions_of_bc(self, regions: 'Regions') -> Dict[int, List[int]]:
        """Returns a dict with the regions for every boundary condition.

        Parameters
        ----------
        regions : Regions
            A regions object.

        Returns
        -------
        regions_oc_bc : Dict[int,List[int]]
            Dict with the key being the ID of a boundary condition and the value being e list of IDs of the regions that
            have this boundary condition.
        """
        regions_of_bc = {k: [] for k in self.get_ids()}
        for key in regions.get_keys():
            bc_id = regions.get_regi(key).bc
            if bc_id is not None:
                regions_of_bc[bc_id].append(key)
        return regions_of_bc

    def shrink(self, msh: 'Mesh', sft: 'ShapeFunction', regions: 'Regions', matrix: coo_matrix, rhs: coo_matrix,
               integration_order: int = 1) -> Tuple[coo_matrix, coo_matrix]:
        """
        Shrinks the system of equations.

        Parameters
        ----------
        msh : Mesh
            The mesh object
        sft : ShapeFunction
            A shape function object
        regions : Regions
            All regions
        matrix : coo_matrix
            Unconstrained FE matrix
        rhs : coo_matrix
            Unconstrained right hand side vector
        integration_order : int, optional
            The order of numerical integration, if needed. Default is 1

        Returns
        -------
        matrix_shrink : coo_matrix
            Constrained FE matrix
        g_shrunk : coo_matrix
            Constrained right hand side vector

        Notes
        -----
        The shrunk system can be solved. To get the final solution, the method inflate has to be used!

        See Also
        --------
        source.bdrycond.BdryCond.BdryCond.inflate : The method to get the final solution.

        Examples
        --------
        Suppose you have msh, regions, and bc. The former two one can get from Geometry. The latter one can contain any
        boundary conditions.

        >>> sft = TriCartesianNodalShapeFunction(msh)
        >>> matrix = sft.divgrad_operator(1)
        >>> q = sft.load_vector(0)
        >>> matrix_shrink, q_shrink = bc.shrink(msh, sft, regions, matrix, rhs)
        >>> u_shrink = spsolve(matrix_shrink, q_shrink)  # From scipy.sparse.linalg
        >>> u = bc.inflate(msh, sft, regions, u_shrink)  # Satisfies matrix * u = rhs
        """
        if not isinstance(matrix, coo_matrix):
            if isinstance(matrix, (csc_matrix, csr_matrix, lil_matrix)):
                matrix = matrix.tocoo()
            else:
                raise ValueError("Matrix is not a coo_matrix")
        if not isinstance(rhs, coo_matrix):
            if isinstance(rhs, (csc_matrix, csr_matrix, lil_matrix)):
                rhs = rhs.tocoo()
            else:
                raise ValueError("Right-hand-side is not a coo_matrix")

        if matrix.shape[0] != rhs.shape[0]:
            raise Exception("The dimensions of matrix and right-hand-side do not match.")

        @dataclass
        class Problem:
            """Temporary problem class."""

            mesh: 'Mesh'
            regions: 'Regions'
            boundary_conditions: 'BdryCond'

        matrix_shrink, rhs_shrink, _, _, _ = sft.shrink(matrix, rhs, Problem(msh, regions, self),
                                                        integration_order=integration_order)

        return matrix_shrink, rhs_shrink

    # noinspection PyUnboundLocalVariable
    def inflate(self, msh: 'Mesh', sft: 'ShapeFunction', regions: 'Regions', u: ndarray) -> ndarray:
        r"""Inflate the system of equations.

        Inserts the information from the boundary conditions into the solution. This method has to be called after the
        system of equation got shrunk and solved. Make sure that the parameters msh, sft, and regions are the same as
        used in shrink.

        Parameters
        ----------
        msh : Mesh
            An instance of the mesh.
        sft : ShapeFunction
            An instance of ShapeFunction.
        regions : Regions
            An instance of Regions
        u : ndarray
            A vector. This is the solution of the shrunk system

        Returns
        -------
        u_inflated : ndarray
            The solution of the original system of equations.

        Notes
        -----
        Let the original system of equations be :math:`\mathbf{A}\mathbf{u} = \mathbf{b}` with additional boundary
        conditions. From the shrink method you get the matrix :math:`\mathbf{A}_{\mathrm{shrink}}` and the vector
        :math:`\mathbf{b}_{\mathrm{shrink}}`. Let :math:`\mathbf{u}_{\mathrm{shrink}}` be the solution of this system.
        Use the inflate method to get :math:`\mathrm{u}`.

        See Also
        --------
        source.bdrycond.BdryCond.BdryCond.shrink : The method to shrink the system.

        Examples
        --------
        Suppose you have msh, regions, and bc. The former two one can get from Geometry. The latter one can contain any
        boundary conditions.

        >>> sft = TriCartesianNodalShapeFunction(msh)
        >>> matrix = sft.divgrad_operator(1)
        >>> q = sft.load_vector(0)
        >>> matrix_shrink, q_shrink = bc.shrink(msh, sft, regions, matrix, rhs)
        >>> u_shrink = spsolve(matrix_shrink, q_shrink)  # From scipy.sparse.linalg
        >>> u = bc.inflate(msh, sft, regions, u_shrink)  # Satisfies matrix * u = rhs
        """

        @dataclass
        class Problem:
            """Temporary problem class."""

            mesh: 'Mesh'
            regions: 'Regions'
            boundary_conditions: 'BdryCond'

        return sft.inflate(u, Problem(msh, regions, self))
