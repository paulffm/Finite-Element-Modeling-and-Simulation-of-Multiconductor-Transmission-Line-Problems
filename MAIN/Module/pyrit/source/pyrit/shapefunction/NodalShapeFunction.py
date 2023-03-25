#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contains the abstract class NodalShapeFunction.

.. sectionauthor:: Christ, Bundschuh
"""

from abc import abstractmethod
from typing import Union, Callable, Tuple, TYPE_CHECKING, List, Dict, Any
import logging

import numpy as np
from numpy import ndarray
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack, bmat

from . import ShapeFunction

if TYPE_CHECKING:
    from pyrit.material import MatProperty, Materials
    from pyrit.region import Regions
    from pyrit.bdrycond import BdryCond, BCDirichlet
    from pyrit.excitation import Excitations
    from pyrit.mesh import Mesh
    from .ShapeFunction import ShrinkInflateProblemShapeFunction


class NodalShapeFunction(ShapeFunction):
    """Abstract class for nodal shape functions.

    See Also
    --------
    pyrit.shapefunction.ShapeFunction : Abstract super class
    pyrit.mesh : Abstract underlying Mesh object on whose entities the SFT parameters are computed.
    """

    @abstractmethod
    def __init__(self, mesh, dim: int, allocation_size: int):
        """
        Constructor for the abstract class NodalShapeFunction.

        Parameters
        ----------
        mesh : Mesh
            Representing the underlying mesh on which the NodalShapeFunction is
            defined.
        dim : int
            Dimensionality of the shape function.
        allocation_size : int
            Size used for vector allocation in matrix creation.

        Returns
        -------
        None
        """
        super().__init__(mesh, dim, allocation_size)

    @abstractmethod
    def divgrad_operator(self, *material: Union[Callable[..., float],
                                                ndarray, float, Tuple['Regions', 'Materials', 'MatProperty']],
                         integration_order: int = 1) -> coo_matrix:
        r"""
        Compute the discrete version of the div-grad operator as a matrix.

        Returns the div-grad-matrix, which is the discrete representation of
        the weak formulation

        .. math::

            (\mathrm{material} \cdot \mathrm{grad}\,u, \mathrm{grad}\,u')

        for the given nodal shape function object.
        Use e.g. to solve an electrostatic problem for the potential solution on the nodes.

        Parameters
        ----------
        material : Union[Callable[..., float], ndarray, float, Tuple[Regions, Materials, MatProperty]]
            Information of the material parameter inside the domain. Needs to be strictly positive > 0 This can be:

                - A function that gives a value for every point inside the domain with as many input arguments as the \
                    dimension of the domain, e.g. f(x,y,z)
                - A (T,) array that has one value for every element (T elements)
                - A single value if the material has everywhere the same value
                - A tuple (Regions, Materials, MatProperty). MatProperty has to be the class of the property that \
                    shall be considered, e.g. Permittivity

        integration_order: int
            If a function is given or the material of one region depends on space, a numerical integration is used to
            calculate the values in the matrix. This is the order of the numerical integration.

        Returns
        -------
        div_grad : (N,N) coo_matrix
            Sparse matrix representing the discrete div-grad operator for nodal shape functions on the underlying mesh.
            (N is the number of nodes in the mesh)

        See Also
        --------
        pyrit.geometry.Geometry: To get the Regions and Materials object
        scipy.sparse: For the coo_matrix that is returned
        """

    @abstractmethod
    def mass_matrix(self, *material: Union[Callable[..., float],
                                           ndarray, float, Tuple['Regions', 'Materials', 'MatProperty']],
                    integration_order: int = 1) -> coo_matrix:
        """
        Compute the discrete version of a material multiplication as a matrix.

        Returns the nodal mass-matrix, which is the discrete representation of the weak formulation

        .. math::     (material * u, u')

        for the given nodal shape function object.

        Parameters
        ----------
        material : Union[Callable[..., float], ndarray, float, Tuple[Regions, Materials, MatProperty]]
            Information of the material parameter inside the domain. Needs to be strictly positive > 0. This can be:

                - A function that returns a scalar value for every point inside the domain with as many input \
                    arguments as the dimension of the domain, e.g. f(x,y,z), f(x,y) or f(x)
                - A (T,) array that has one value for every element (T elements)
                - A single value if the material has everywhere the same value
                - A tuple (Regions, Materials, MatProperty). MatProperty has to be the class of the property that \
                    shall be considered, e.g. Permittivity

        integration_order: int
            If a function is given or the material of one region depends on space, a numerical integration is used to
            calculate the values in the matrix. This is the order of the numerical integration.

        Returns
        -------
        node_mass : (N,N) coo_matrix
            Sparse matrix for nodal shape functions on the underlying mesh. (N is the number of nodes in the mesh)

        See Also
        --------
        pyrit.geometry.Geometry: To get the Regions and Materials object
        scipy.sparse: For the coo_matrix that is returned
        """

    @abstractmethod
    def load_vector(self, *load: Union[Callable[..., float], ndarray, float, Tuple['Regions', 'Excitations']],
                    integration_order: int = 1) -> coo_matrix:
        """Compute the discrete version of a node-wise excitation as a vector.

        Returns the nodal load vector, typically the equivalent to a
        right-hand-side excitation in a PDE.
        Use e.g. to calculate the discrete version of a charge distribution in
        an electrostatic problem with the potential solution on the nodes.

        Parameters
        ----------
        load : Union[Callable[..., float], ndarray, float, Tuple[Regions, Excitations]]
            Information of the excitations in the domain. This can be:

                - A function that gives a value for every point inside the domain with as many input arguments as the \
                    dimension of the domain, e.g. f(x,y,z)
                - A (T,) array that has one value for every element (T elements)
                - A single value if the excitation has everywhere the same value
                - A tuple (Regions, Excitations)

        integration_order: int
            If a function is given or the excitation of one region depends on space, a numerical integration is used to
            calculate the values in the vector. This is the order of the numerical integration.

        Returns
        -------
        node_load : (N,1) coo_matrix
            Nodal load vector. (N is the number of nodes in the mesh)

        Raises
        ------
        ValueError :
            if input load does not match one of the above formats or input integration_order is not implemented or
            has an invalid format.

        """

    @abstractmethod
    def gradient(self, val2node: ndarray) -> ndarray:
        r"""Compute the gradient of the given scalar field.

        Returns the result of

        .. math::

            \mathrm{grad}(\text{val2node})

        which means it translates a scalar function on the nodes to a vector field on the elements.
        Use e.g. to calculate the electric field from a given potential during postprocessing.

        Parameters
        ----------
        val2node : (N,) ndarray
            Nodal values of a scalar function.

        Returns
        -------
        grad : {(T,), (T,2), (T,3)} ndarray
            Matrix assigning a vector to each element, dimension depending on
            dimension of mesh (1D, 2D, 3D)

        """

    @abstractmethod
    def neumann_term(self, *args: Union[Tuple['Regions', 'BdryCond'],
                                        Tuple[ndarray, Union[Callable[..., float], ndarray, float]]],
                     integration_order: int = 1) -> coo_matrix:
        r"""Compute the Neumann term on a part of the boundary (see the notes).

        Parameters
        ----------
        args : Union[Tuple[Regions, BdryCond], Tuple[ndarray, Union[Callable[..., float], ndarray, float]]]
            The following possibilities are:

            - **Regions, BdyCond**: Will search for instances of BCNeumann in BdryCond and use the value of these.
            - **ndarray, Union[Callable[..., float], ndarray, float**: Indices of the boundary elements and the value of
              :math:`g` at the elements. This can be a function, an array with one value per boundary element or a
              constant.

        integration_order : int, optional
            The integration order for the case that a numerical integration is necessary. Default is 1

        Returns
        -------
        q : coo_matrix
            The vector :math:`\mathbf{q}` for the Neumann term.

        Notes
        -----
        The Neumann term for nodal basis functions is a vector :math:`\mathbf{q}` where the elements are defined as

        .. math::

            q_i = \int_{\partial\Omega}g N_i\;\mathrm{d} S\,.
        """

    @abstractmethod
    def robin_terms(self, regions: 'Regions', boundary_condition: 'BdryCond', integration_order: int = 1) -> \
            Tuple[coo_matrix, coo_matrix]:
        r"""Compute the matrix :math:`\mathbf{S}` and vector :math:`\mathbf{q}` arising from Robin boundary conditions

        See the Notes section for detailed information about the matrix and vector.

        Parameters
        ----------
        regions: Regions
            An instance of Regions
        boundary_condition: BdryCond
            An instance of BdryCond
        integration_order : int, optional
            The integration order for the case that a numerical integration is necessary. Default is 1

        Returns
        -------
        S : coo_matrix
            The matrix for the left hand side. A sparse matrix of dimension (E,E), with E being the number of elements.
        q : coo_matrix
            The vector for the right hand side. A sparse vector of dimension (E,1), with E being the number of elements.

        Notes
        -----
            For nodal basis functions and the Robin boundary condition given as
            :math:`\alpha\phi + \beta \frac{\partial\phi}{\partial n} = g`, the elements of the matrix
            :math:`\mathbf{S}` and vector :math:`\mathbf{q}` are defined as

            .. math::

                S_{ij} &= \int_{\partial\Omega}\frac\alpha\beta N_j N_i\;\mathrm{d} S\,,\\
                q_i &= \int_{\partial\Omega}\frac1\beta gN_i\;\mathrm{d} S\,.
        """

    def shrink(self, matrix: coo_matrix, rhs: coo_matrix, problem: 'ShrinkInflateProblemShapeFunction',
               integration_order: int = 1) -> Tuple['coo_matrix', 'coo_matrix', ndarray, int, Dict['str', Any]]:
        # Compute needed data
        bc_by_type = problem.boundary_conditions.dict_of_boundary_condition()
        regions_of_bc = problem.boundary_conditions.regions_of_bc(problem.regions)

        # Initialize some variables
        num_lagrange_mul_binary = 0  # Number of lagrange multipliers for binary boundary conditions
        num_lagrange_mul_floating = 0  # Number fo lagrange multipliers for floating boundary conditions
        ind_on_dir, ind_not_on_dir, val_on_dir = None, None, None

        # For every type of boundary condition, if it exists, shrink for this type

        if bc_by_type['neumann']:
            rhs = self.shrink_neumann(rhs, problem.boundary_conditions, problem.regions, integration_order)

        if bc_by_type['robin']:
            matrix, rhs = self.shrink_robin(matrix, rhs, problem.boundary_conditions, problem.regions,
                                            integration_order)

        if bc_by_type['binary']:
            matrix, rhs, num_lagrange_mul_binary = self.shrink_binary(matrix, rhs, problem.boundary_conditions,
                                                                      bc_by_type['binary'])

        if bc_by_type['floating']:
            matrix, rhs, num_lagrange_mul_floating = self.shrink_floating(matrix, rhs, bc_by_type['floating'],
                                                                          problem.regions, problem.mesh,
                                                                          regions_of_bc)

        if bc_by_type['dirichlet']:
            out = self.shrink_dirichlet(matrix, rhs, problem.boundary_conditions, bc_by_type['dirichlet'],
                                        problem.regions, problem.mesh, regions_of_bc)
            matrix, rhs, ind_on_dir, ind_not_on_dir, val_on_dir = out

        # Prepare output data
        num_lagrange_mul = num_lagrange_mul_binary + num_lagrange_mul_floating
        support_data = {'num_lagrange_mul_binary': num_lagrange_mul_binary,
                        'num_lagrange_mul_floating': num_lagrange_mul_floating,
                        'indices_on_dirichlet': ind_on_dir,
                        'indices_not_on_dirichlet': ind_not_on_dir,
                        'values_on_dirichlet': val_on_dir}

        indices_not_dof = np.array([], dtype=int)

        if ind_on_dir is not None:
            indices_not_dof = ind_on_dir

        return matrix, rhs, indices_not_dof, num_lagrange_mul, support_data

    def inflate(self, solution: ndarray, problem: 'ShrinkInflateProblemShapeFunction',
                support_data: Dict[str, Any] = None) -> ndarray:
        # Compute needed data
        bc_by_type = problem.boundary_conditions.dict_of_boundary_condition()

        # Check if support data is given and complete
        if support_data:
            required_keys = ['num_lagrange_mul_binary', 'num_lagrange_mul_floating',
                             'indices_on_dirichlet', 'indices_not_on_dirichlet', 'values_on_dirichlet']
            for key in required_keys:
                if key not in support_data.keys():
                    # pylint: disable=logging-too-many-args
                    logging.warning(
                        "The given support data is incomplete. The key '%s' is missing. It will be computed.", key)
                    support_data = self.compute_support_data(bc_by_type, len(solution), problem)
                    break
        else:
            support_data = self.compute_support_data(bc_by_type, len(solution), problem)

        # Read the support data
        num_lagrange_mul_binary = support_data['num_lagrange_mul_binary']
        num_lagrange_mul_floating = support_data['num_lagrange_mul_floating']
        ind_on_dir = support_data['indices_on_dirichlet']
        ind_not_on_dir = support_data['indices_not_on_dirichlet']
        val_on_dir = support_data['values_on_dirichlet']

        # For every type of boundary condition, if it exists, inflate this type

        if bc_by_type['dirichlet']:
            solution_inflated = self.inflate_dirichlet(solution, ind_on_dir, val_on_dir, ind_not_on_dir)
        else:
            solution_inflated = solution.copy()

        if bc_by_type['floating']:
            solution_inflated = self.inflate_floating(solution_inflated, num_lagrange_mul_floating)

        if bc_by_type['binary']:
            solution_inflated = self.inflate_binary(solution_inflated, num_lagrange_mul_binary)

        return solution_inflated

    # region shrink methods

    def shrink_neumann(self, rhs: coo_matrix, boundary_conditions: 'BdryCond', regions: 'Regions',
                       integration_order: int = 1) -> coo_matrix:
        """Shrink with Neumann boundary conditions.

        Parameters
        ----------
        rhs : coo_matrix
            The right-hand-side of the system of equations
        boundary_conditions : BdryCond
            A boundary condition object
        regions : Regions
            A regions object
        integration_order : int, optional
            The integration order for subsequent methods. Default is 1

        Returns
        -------
        rhs_shrink : coo_matrix
            The shrunk right-hand-side.
        """
        rhs_shrink = rhs + self.neumann_term(regions, boundary_conditions, integration_order=integration_order)
        return rhs_shrink.tocoo()

    def shrink_robin(self, matrix: coo_matrix, rhs: coo_matrix, boundary_conditions: 'BdryCond', regions: 'Regions',
                     integration_order: int = 1) -> Tuple[coo_matrix, coo_matrix]:
        """Shrink with Robin boundary conditions.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the system of equations.
        rhs : coo_matrix
            The right-hand-side of the system of equations.
        boundary_conditions : BdryCond
            A boundary condition object.
        regions : Regions
            A regions object
        integration_order : int, optional
            The integration order for subsequent methods. Default is 1.

        Returns
        -------
        matrix_shrink : coo_matrix
            The shrunk matrix
        rhs_shrunk : coo_matrix
            The shrunk right-hand-side
        """
        matrix_neumann, vector_neumann = self.robin_terms(regions, boundary_conditions,
                                                          integration_order=integration_order)
        matrix_shrink = matrix + matrix_neumann
        rhs_shrink = rhs + vector_neumann
        return matrix_shrink.tocoo(), rhs_shrink.tocoo()

    def shrink_binary(self, matrix: coo_matrix, rhs: coo_matrix, boundary_conditions: 'BdryCond',
                      binary_boundary_conditions: List[int]):
        """Shrink with binary boundary conditions.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the system of equations.
        rhs : coo_matrix
            The right-hand-side of the system of equations.
        boundary_conditions : BdryCond
            A boundary condition object.
        binary_boundary_conditions : List[int]
            A list of the IDs of binary boundary conditions.

        Returns
        -------
        matrix_shrink : coo_matrix
            The shrunk matrix
        rhs_shrunk : coo_matrix
            The shrunk right-hand-side
        num_lagrange_mul : int
            The number of lagrange multipliers
        """
        lines, columns, values = self._calc_values_for_binary(boundary_conditions, binary_boundary_conditions)

        # Build the matrix B, the zeros vector for the right-hand-side and the zero block matrix
        binary_matrix = coo_matrix(
            (np.concatenate(values), (np.concatenate(lines), np.concatenate(columns))))
        binary_matrix.resize((binary_matrix.shape[0], matrix.shape[1]))
        num_lagrange_mul = binary_matrix.shape[0]  # Number of lagrange multiplicators
        binary_vector = coo_matrix((num_lagrange_mul, 1))
        bottom_right_matrix = coo_matrix((num_lagrange_mul, num_lagrange_mul))

        # Update system matrix and rhs vector
        matrix = hstack((vstack((matrix, binary_matrix)),
                         vstack((binary_matrix.transpose(), bottom_right_matrix))), format='coo')
        rhs = vstack((rhs, binary_vector), format='coo')
        return matrix, rhs, num_lagrange_mul


    def shrink_floating(self, matrix: coo_matrix, rhs: coo_matrix, floating_boundary_conditions: List[int],
                        regions: 'Regions', mesh: 'Mesh', regions_of_bc: Dict[int, List[int]]) \
            -> Tuple[coo_matrix, coo_matrix, List[ndarray], int, ndarray]:
        """Shrink with floating boundary conditions.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the system of equations.
        rhs : coo_matrix
            The right-hand-side of the system of equations.
        floating_boundary_conditions : List[int]
            A list of the IDs of floating boundary conditions.
        regions : Regions
            A regions object
        mesh : Mesh
            A mesh object
        regions_of_bc : Dict[int, List[int]]
            A dictionary with the key being the ID of a boundary conditions and the value being a list of the IDs of
            regions having this boundary condition. See :py:meth:`.BdryCond.regions_of_bc`.

        Returns
        -------
        matrix_shrink : coo_matrix
            The shrunk matrix
        rhs_shrink : coo_matrix
            The shrunk right-hand-side
        num_lagrange_mul : int
            Number of lagrange multipliers
        """
        indices_all_nodes = self._calc_values_for_floating(floating_boundary_conditions, mesh,
                                                           regions, regions_of_bc)

        b_matrices = []  # List of diagonal matrices
        for nodes in indices_all_nodes:
            nodes: np.ndarray
            size = nodes.size - 1
            matrix_slave = coo_matrix((np.ones(size), (np.arange(size), nodes[:-1])), shape=(size, matrix.shape[0]))
            matrix_master = coo_matrix((np.ones(size), (np.arange(size), int(nodes[-1]) * np.ones(size, dtype=int))),
                                       shape=(size, matrix.shape[0]))
            b_matrices.append(matrix_slave - matrix_master)

        b_matrix: coo_matrix = vstack(b_matrices, format='coo')

        num_lagrange_mul = b_matrix.shape[0]

        shrink_matrix = bmat([[matrix, b_matrix.transpose()],
                              [b_matrix, None]], format='coo')
        shrink_rhs = vstack([rhs, coo_matrix((num_lagrange_mul, 1))])

        return shrink_matrix, shrink_rhs, num_lagrange_mul

    def shrink_dirichlet(self, matrix: coo_matrix, rhs: coo_matrix, boundary_conditions: 'BdryCond',
                         dirichlet_boundary_conditions: List['BCDirichlet'], regions: 'Regions', mesh: 'Mesh',
                         regions_of_bc: Dict[int, List[int]]) \
            -> Tuple[coo_matrix, coo_matrix, ndarray, ndarray, ndarray, ndarray]:
        """Shrink with dirichlet boundary conditions.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the system of equations.
        rhs : coo_matrix
            The right-hand-side of the system of equations.
        boundary_conditions : BdryCond
            A boundary condition object
        dirichlet_boundary_conditions : List[BCDirichlet]
            A list of dirichlet boundary conditions.
        regions : Regions
            A regions object.
        mesh : Mesh
            A mesh object.
        regions_of_bc : Dict[int, List[int]]
            A dictionary with the key being the ID of a boundary conditions and the value being a list of the IDs of
            regions having this boundary condition. See :py:meth:`.BdryCond.regions_of_bc`.

        Returns
        -------
        matrix_shrink : coo_matrix
            The shrunk matrix
        rhs_shrink : coo_matrix
            The shrunk right-hand-side
        indices_on_dirichlet : ndarray
            An array containing the indices of nodes that are on a dirichlet boundary condition. These indices are valid
            for the input matrix.
        indices_not_on_dirichlet : ndarray
            An array containing the indices of nodes that are not on a dirichlet boundary condition. These indices are
            valid for the input matrix.
        values_on_dirichlet : ndarray
            An array with the value of every node on a dirichlet boundary condition.
        """
        matrix = matrix.tocsr()
        rhs = rhs.tocsr()

        indices_all_nodes, values_all_nodes = self._calc_values_for_dirichlet(boundary_conditions,
                                                                              dirichlet_boundary_conditions, mesh,
                                                                              regions, regions_of_bc)
        values_on_dirichlet = np.concatenate(values_all_nodes)
        indices_on_dirichlet = np.concatenate(indices_all_nodes)

        indices_not_on_dirichlet = np.setdiff1d(np.arange(matrix.shape[0]), indices_on_dirichlet)

        # Update system matrix and right hand side vector
        rhs = rhs[indices_not_on_dirichlet] - matrix[
            np.ix_(indices_not_on_dirichlet, indices_on_dirichlet)] * csr_matrix(values_on_dirichlet)
        matrix = matrix[np.ix_(indices_not_on_dirichlet, indices_not_on_dirichlet)]

        return matrix.tocoo(), rhs.tocoo(), indices_on_dirichlet, indices_not_on_dirichlet, values_on_dirichlet.reshape(
            -1)

    # endregion

    # region inflate methods

    @staticmethod
    def inflate_dirichlet(solution: ndarray, indices_on_dirichlet: ndarray, values_on_dirichlet: ndarray,
                          indices_not_on_dirichlet: ndarray) -> ndarray:
        """Inflate Dirichlet boundary conditions.

        Parameters
        ----------
        solution : ndarray
            The solution of the shrunk system.
        indices_on_dirichlet : ndarray
            An array containing the indices of nodes that are on a dirichlet boundary condition.
            See :py:meth:`shrink_dirichlet`
        values_on_dirichlet : ndarray
            An array with the value of every node on a dirichlet boundary condition. See :py:meth:`shrink_dirichlet`
        indices_not_on_dirichlet : ndarray
            An array containing the indices of nodes that are not on a dirichlet boundary condition.
            See :py:meth:`shrink_dirichlet`

        Returns
        -------
        solution_inflated : ndarray
            The inflated solution.
        """
        solution_inflated = np.empty(len(indices_on_dirichlet) + len(indices_not_on_dirichlet), dtype=solution.dtype)
        solution_inflated[indices_not_on_dirichlet] = solution
        solution_inflated[indices_on_dirichlet] = values_on_dirichlet

        return solution_inflated

    @staticmethod
    def inflate_floating(solution: ndarray, num_lagrange_mul: int) -> ndarray:
        """Inflate floating boundary conditions.

        Parameters
        ----------
        solution : ndarray
            The solution of the shrunk system.
        num_lagrange_mul : int
            The number of lagrange multipliers. See :py:meth:`shrink_floating`

        Returns
        -------
        solution_inflated : ndarray
            The inflated solution.
        """
        solution_inflated = solution[:len(solution) - num_lagrange_mul]
        return solution_inflated

    @staticmethod
    def inflate_binary(solution: ndarray, num_lagrange_mul: int) -> ndarray:
        """Inflate binary boundary conditions.

        Parameters
        ----------
        solution : ndarray
            The solution of the shrunk system.
        num_lagrange_mul : int
            The number of lagrange multipliers. See :py:meth:`shrink_binary`

        Returns
        -------
        solution_inflated : ndarray
            The inflated solution.
        """
        solution_inflated = solution[:len(solution) - num_lagrange_mul]
        return solution_inflated

    # endregion

    # region helper methods

    # noinspection PyUnresolvedReferences
    @staticmethod
    def _calc_values_for_binary(boundary_conditions: 'BdryCond', keys: List[int]) \
            -> Tuple[List[int], List[int], List[float]]:
        """Caluclate some data needed for shrinking binary boundary conditions.

        Parameters
        ----------
        boundary_conditions : BdryCond
            A boundary condition object
        keys : List[int]
            A list of the IDs of binary boundary conditions.

        Returns
        -------
        lines : List[int]
            Line numbers
        columns : List[int]
            Column numbers
        values : List[float]
            Values of the binary boundary conditions.
        """
        lines = []
        columns = []
        values = []
        lines_previous = 0
        for key in keys:
            bc_binary = boundary_conditions.get_bc(key)

            lines.append(np.arange(len(bc_binary.primary_nodes)) + lines_previous)
            columns.append(bc_binary.primary_nodes)
            values.append(-bc_binary.value * np.ones(len(bc_binary.primary_nodes)))

            lines.append(np.arange(len(bc_binary.replica_nodes)) + lines_previous)
            columns.append(bc_binary.replica_nodes)
            values.append(np.ones(len(bc_binary.replica_nodes)))
            lines_previous += len(bc_binary.primary_nodes)
        return lines, columns, values

    # noinspection PyUnresolvedReferences
    @staticmethod
    def _calc_values_for_floating(keys: List[int], msh: 'Mesh', regions: 'Regions',
                                  regions_of_bc: Dict[int, List[int]]) -> List[ndarray]:
        """Compute some values needed for shrinking floating boundary conditions.

        Parameters
        ----------
        keys : List[int]
            A list of the indices of floating boundary conditions.
        msh : Mesh
            A mesh object
        regions : Regions
            A regions object
        regions_of_bc : Dict[int, List[int]]
            A dictionary with the key being the ID of a boundary conditions and the value being a list of the IDs of
            regions having this boundary condition. See :py:meth:`.BdryCond.regions_of_bc`.

        Returns
        -------
        indices_all_nodes : List[ndarray]
            A list of arrays. Each array contains the nodes on one floating boundary condition. Especially, there are
            the nodes of different regions that have the same boundary condition.
        """
        indices_all_nodes = []  # Indices of all nodes on a floating boundary condition
        # Get one array of the node indices for every floating boundary condition
        for key in keys:
            regions_of_floating = regions_of_bc[key]

            # Get all node indices for this boundary condition
            indices_nodes_list = []
            for regi in regions_of_floating:
                dimension_of_region = regions.get_regi(regi).dim
                if dimension_of_region == 1:
                    indices_entities = np.where(msh.edge2regi == regi)
                    indices_nodes_list.append(np.unique(msh.edge2node[indices_entities, :]))
                elif dimension_of_region == 2:
                    indices_entities = np.where(msh.edge2regi == regi)
                    indices_nodes_list.append(np.unique(msh.edge2node[indices_entities, :]))
            indices_all_nodes.append(np.unique(np.concatenate(indices_nodes_list)))

        return indices_all_nodes

    # noinspection PyUnresolvedReferences
    @staticmethod
    def _calc_values_for_dirichlet(boundary_conditions: 'BdryCond', keys: List[int], msh: 'Mesh', regions: 'Regions',
                                   regions_of_bc: Dict[int, List[int]]) -> Tuple[List[int], List[float]]:
        """Calculate some data needed for shrinking and inflating with Dirichlet boundary conditions.

        Parameters
        ----------
        boundary_conditions : BdryCond
            A boundary condition object.
        keys : List[int]
            The IDs of the Dirichlet boundary conditions.
        msh : Mesh
            A mesh object
        regions : Regions
            A regions object.
        regions_of_bc : Dict[int, List[int]]
            A dictionary with the key being the ID of a boundary conditions and the value being a list of the IDs of
            regions having this boundary condition. See :py:meth:`.BdryCond.regions_of_bc`.

        Returns
        -------
        indices_all_nodes : List[int]
            A list with indices.
        values_all_nodes : List[float]
            A list with values.
        """
        indices_all_nodes = []
        values_all_nodes = []
        for key in keys:
            bc_dir = boundary_conditions.get_bc(key)
            regions_of_dirichlet = regions_of_bc[key]

            # Get all node indices for this boundary condition
            indices_nodes_list = []
            for regi in regions_of_dirichlet:
                dimension_of_region = regions.get_regi(regi).dim
                if dimension_of_region == 0:
                    indices_nodes_list.append(np.where(msh.node2regi == regi)[0])
                if dimension_of_region == 1:
                    indices_entities = np.where(msh.edge2regi == regi)[0]
                    indices_nodes_list.append(np.unique(msh.edge2node[indices_entities, :]))
                elif dimension_of_region == 2:
                    indices_entities = np.where(msh.elem2regi == regi)[0]
                    indices_nodes_list.append(np.unique(msh.elem2node[indices_entities, :]))
            indices_nodes = np.unique(np.concatenate(indices_nodes_list))
            indices_all_nodes.append(indices_nodes)

            if bc_dir.is_constant:
                if isinstance(bc_dir.value, np.ndarray):
                    if bc_dir.value.size != len(indices_nodes):
                        raise Exception("Dimensions dont match")
                    values_all_nodes.append(bc_dir.value.reshape(-1, 1))
                elif isinstance(bc_dir.value, (float, int, complex)):
                    values_all_nodes.append(bc_dir.value * np.ones((len(indices_nodes), 1)))
                else:
                    raise ValueError("Type not supported")
            elif bc_dir.is_time_dependent:
                values_all_nodes.append(bc_dir.value() * np.ones((len(indices_nodes), 1)))
            else:
                raise NotImplementedError

        # When Dirichlet boundary conditions intersect each other
        indices_tmp = indices_all_nodes[0]
        for k in range(1, len(indices_all_nodes)):
            intersection = np.intersect1d(indices_tmp, indices_all_nodes[k])
            if len(intersection) > 0:
                indices_cut = []
                for m in intersection:
                    indices_cut.append(np.where(indices_all_nodes[k] == m)[0][0])
                indices_keep = np.setdiff1d(np.arange(len(indices_all_nodes[k])), indices_cut)
                indices_all_nodes[k] = indices_all_nodes[k][indices_keep]
                values_all_nodes[k] = values_all_nodes[k][indices_keep]
            indices_tmp = np.concatenate(indices_all_nodes[:k + 1])

        return indices_all_nodes, values_all_nodes

    def compute_support_data(self, dict_of_bc: Dict[int, List[int]], size_solution: int,
                             problem: 'ShrinkInflateProblemShapeFunction') -> Dict[str, Any]:
        """Compute support data needed for inflating.

        Parameters
        ----------
        dict_of_bc : Dict[int, List[int]]
            A dictionary with the key being the ID of a boundary conditions and the value being a list of the IDs of
            regions having this boundary condition. See :py:meth:`.BdryCond.regions_of_bc`.
        size_solution : int
            Size of the solution vector.
        problem : ShrinkInflateProblemShapeFunction
            A problem object.

        Returns
        -------
        support_data : Dict[str, Any]
            Support data needed for inflating.
        """
        regions_of_bc = problem.boundary_conditions.regions_of_bc(problem.regions)
        num_lagrange_mul_binary = 0
        num_lagrange_mul_floating = 0
        indices_on_dirichlet, indices_not_on_dirichlet, values_on_dirichlet = None, None, None

        if dict_of_bc['binary']:
            lines, _, _ = self._calc_values_for_binary(problem.boundary_conditions, dict_of_bc['binary'])
            num_lagrange_mul_binary = np.max(lines) + 1

        if dict_of_bc['floating']:
            indices_all_nodes = self._calc_values_for_floating(dict_of_bc['floating'], problem.mesh, problem.regions,
                                                               regions_of_bc)
            num_lagrange_mul_floating = sum(nodes.size - 1 for nodes in indices_all_nodes)

        if dict_of_bc['dirichlet']:
            indices_all_nodes_dirichlet, values_all_nodes = self._calc_values_for_dirichlet(problem.boundary_conditions,
                                                                                            dict_of_bc['dirichlet'],
                                                                                            problem.mesh,
                                                                                            problem.regions,
                                                                                            regions_of_bc)
            values_on_dirichlet = np.concatenate(values_all_nodes).reshape(-1)
            indices_on_dirichlet = np.concatenate(indices_all_nodes_dirichlet)

            indices_not_on_dirichlet = np.setdiff1d(np.arange(size_solution + len(indices_on_dirichlet)),
                                                    indices_on_dirichlet)

        support_data = {'num_lagrange_mul_binary': num_lagrange_mul_binary,
                        'num_lagrange_mul_floating': num_lagrange_mul_floating,
                        'indices_on_dirichlet': indices_on_dirichlet,
                        'indices_not_on_dirichlet': indices_not_on_dirichlet,
                        'values_on_dirichlet': values_on_dirichlet}

        return support_data

    # endregion
