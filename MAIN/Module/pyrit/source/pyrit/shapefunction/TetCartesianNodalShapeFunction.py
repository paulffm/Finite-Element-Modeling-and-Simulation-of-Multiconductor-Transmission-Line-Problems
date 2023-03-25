# coding=utf-8
"""
Created on 21 June 2021

File containing the class TetCartesianNodalShapeFunction

.. sectionauthor:: christ
"""
from inspect import signature
from typing import Union, Callable, Tuple  # , TYPE_CHECKING
from warnings import warn

import numpy as np
from scipy.sparse import coo_matrix
from numba import prange

from pyrit.bdrycond import BdryCond, BCNeumann, BCRobin

from pyrit.mesh import TetMesh
from pyrit.material import MatProperty
from pyrit.region import Regions
from pyrit.material import Materials
from pyrit.excitation import Excitations
from pyrit.toolbox import QuadratureToolbox

from . import NodalShapeFunction


def calc_evaluation_points_element(local_coordinates: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray,
                                   offset: np.ndarray,
                                   elements: np.ndarray) -> np.ndarray:
    r"""
    Compute the coordinates at which a function needs to be evaluated if numerical integration is needed.

    Returns a three-dimensional array with the evaluation coordinates for all elements specified. The returned array is
    of dimensions (N,E,T) with [n,e,0] being the x coordinates and [n,e,1] being the y coordinates. Mathematically
    the transformation can be expressed as :math:`(x,\,y)^T  = \mathbf{B}\,(\hat{x},\,\hat{y})^T + offset`.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements
        T       Number of mesh elements
        ======  =======

    Parameters
    ----------
    local_coordinates : np.ndarray
        (N,3) array. x and y coordinates for each evaluation point.
    a : np.ndarray
        (T,3) array. Difference vector of node coordinates.
    b : np.ndarray
        (T,3) array. Difference vector of node coordinates.
    c : np.ndarray
        (T,3) array. Difference vector of node coordinates.
    offset : np.ndarray
        (T,3) array. Offset transformation to global coordinates.
    elements : np.ndarray
        (E,) array. Array with indices of considered elements.

    Returns
    -------
    coord : np.ndarray
        (N,E,3) array.
    """
    matrix_b = np.empty((len(elements), 3, 3))
    matrix_b[:, :, 0] = a[elements, :]
    matrix_b[:, :, 1] = b[elements, :]
    matrix_b[:, :, 2] = c[elements, :]

    coord = np.empty((len(elements), local_coordinates.shape[0], 3))
    for counter in prange(len(elements)):
        coord[counter, :, :] = (matrix_b[counter] @ local_coordinates.T).T + offset[counter]
    return coord


def eval_inhom_lin_iso(fun: Callable[[float, float, float], float], evaluation_points: np.ndarray) -> np.ndarray:
    """
    Evaluate an inhomogeneous, linear and isotropic function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float, float], float]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,3) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

    Returns
    -------
    evaluations : np.ndarray
        (E,N) array. For every element the N evaluations of the material function at the evaluation points for the
        numerical integration.

    Notes
    -----
    This function has to be called before matrix elements, where a numerical integration of an arbitrary function is
    necessary, can be computed.
    """
    values = np.zeros((evaluation_points.shape[0], evaluation_points.shape[1]))
    for k in range(evaluation_points.shape[0]):
        for kk in range(evaluation_points.shape[1]):
            values[k, kk] = fun(evaluation_points[k, kk, 0], evaluation_points[k, kk, 1], evaluation_points[k, kk, 2])
    return values


def calc_divgrad_constant_scalar(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: float, inv_b_trans: np.ndarray, dphi1: np.ndarray, dphi2: np.ndarray,
                                 dphi3: np.ndarray, dphi4: np.ndarray, tetrahedral_volume_ref: float,
                                 det_b: np.ndarray, mesh_elem2node: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: float
        Material value.
    inv_b_trans : np.ndarray
        (E,3,3) array. Inverse and transpose of B-matrix.
    dphi1 : np.ndarray
        (1,3) array. Derivative of first base function.
    dphi2 : np.ndarray
        (1,3) array. Derivative of second base function.
    dphi3 : np.ndarray
        (1,3) array. Derivative of third base function.
    dphi4 : np.ndarray
        (1,3) array. Derivative of fourth base function.
    tetrahedral_volume_ref : float
        Volume of reference tetrahedral.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.

    Returns
    -------
    None
    """
    calc_divgrad_scalar_per_elem(indices, i, j, v, value * np.ones_like(indices), inv_b_trans, dphi1, dphi2, dphi3,
                                 dphi4, tetrahedral_volume_ref, det_b, mesh_elem2node)


def calc_divgrad_scalar_per_elem(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: np.ndarray, inv_b_trans: np.ndarray, dphi1: np.ndarray, dphi2: np.ndarray,
                                 dphi3: np.ndarray, dphi4: np.ndarray, tetrahedral_volume_ref: float,
                                 det_b: np.ndarray, mesh_elem2node: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: np.ndarray
        (E,) array. Material value.
    inv_b_trans : np.ndarray
        (E,3,3) array. Inverse and transpose of B-matrix.
    dphi1 : np.ndarray
        (1,3) array. Derivative of first base function.
    dphi2 : np.ndarray
        (1,3) array. Derivative of second base function.
    dphi3 : np.ndarray
        (1,3) array. Derivative of third base function.
    dphi4 : np.ndarray
        (1,3) array. Derivative of fourth base function.
    tetrahedral_volume_ref : float
        Volume of reference tetrahedral.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.

    Returns
    -------
    None
    """
    i_local = np.zeros((16 * len(indices)), dtype=np.int_)
    j_local = np.zeros_like(i_local, dtype=np.int_)
    v_local = np.zeros_like(i_local, dtype=np.float_)
    indices_ijv = np.zeros_like(i_local, dtype=np.int_)

    for count in prange(len(indices)):  # iterate over each given mesh element
        idx = indices[count]
        matrix_prod = inv_b_trans[:, :, idx] @ np.array([dphi1, dphi2, dphi3, dphi4])[:, :, 0].T
        dg_elem = tetrahedral_volume_ref * det_b[idx] * value[count] * matrix_prod.T @ matrix_prod
        idx_node = mesh_elem2node[idx, :]
        idx_ijv = np.arange(16 * count, 16 * count + 16)
        i_local[idx_ijv] = np.column_stack((idx_node, idx_node, idx_node, idx_node)).flatten()
        j_local[idx_ijv] = np.hstack((idx_node, idx_node, idx_node, idx_node))
        v_local[idx_ijv] = dg_elem.flatten()
        indices_ijv[idx_ijv] = np.arange(16 * idx, 16 * idx + 16)

    i[indices_ijv] = i_local
    j[indices_ijv] = j_local
    v[indices_ijv] = v_local


def calc_divgrad_function_scalar(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: np.ndarray, inv_b_trans: np.ndarray, dphi1: np.ndarray, dphi2: np.ndarray,
                                 dphi3: np.ndarray, dphi4: np.ndarray, det_b: np.ndarray,
                                 mesh_elem2node: np.ndarray, weights: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points per element
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: np.ndarray
        (E, N). Material value.
    inv_b_trans : np.ndarray
        (E,3,3) array. Inverse and transpose of B-matrix.
    dphi1 : np.ndarray
        (1,3) array. Derivative of first base function.
    dphi2 : np.ndarray
        (1,3) array. Derivative of second base function.
    dphi3 : np.ndarray
        (1,3) array. Derivative of third base function.
    dphi4 : np.ndarray
        (1,3) array. Derivative of fourth base function.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.
    weights : np.ndarray
        (N,) array. Weights for numerical integration.

    Returns
    -------
    None
    """
    i_local = np.zeros((16 * len(indices)), dtype=np.int_)
    j_local = np.zeros_like(i_local, dtype=np.int_)
    v_local = np.zeros_like(i_local, dtype=np.float_)
    indices_ijv = np.zeros_like(i_local, dtype=np.int_)

    for count in prange(len(indices)):  # iterate over each given mesh element
        idx = indices[count]

        matrix_prod = inv_b_trans[:, :, idx] @ np.array([dphi1, dphi2, dphi3, dphi4])[:, :, 0].T
        dg_elem = det_b[idx] * np.sum(value[count] * weights) * matrix_prod.T @ matrix_prod
        idx_node = mesh_elem2node[idx, :]
        idx_ijv = np.arange(16 * count, 16 * count + 16)
        i_local[idx_ijv] = np.column_stack((idx_node, idx_node, idx_node, idx_node)).flatten()
        j_local[idx_ijv] = np.hstack((idx_node, idx_node, idx_node, idx_node))
        v_local[idx_ijv] = dg_elem.flatten()
        indices_ijv[idx_ijv] = np.arange(16 * idx, 16 * idx + 16)

    i[indices_ijv] = i_local
    j[indices_ijv] = j_local
    v[indices_ijv] = v_local


def calc_mass_constant_scalar(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                              value: float, tetrahedral_volume_ref: float,
                              det_b: np.ndarray, mesh_elem2node: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: float
        Material value.
    tetrahedral_volume_ref : float
        Volume of reference tetrahedral.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.

    Returns
    -------
    None
    """
    calc_mass_scalar_per_elem(indices, i, j, v, value * np.ones_like(indices), tetrahedral_volume_ref, det_b,
                              mesh_elem2node)


def calc_mass_scalar_per_elem(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                              value: np.ndarray, tetrahedral_volume_ref: float,
                              det_b: np.ndarray, mesh_elem2node: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: np.ndarray.
        (E,) array. Material value.
    tetrahedral_volume_ref : float
        Volume of reference tetrahedral.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.

    Returns
    -------
    None
    """
    i_local = np.zeros((16 * len(indices)), dtype=np.int_)
    j_local = np.zeros_like(i_local, dtype=np.int_)
    v_local = np.zeros_like(i_local, dtype=np.float_)
    indices_ijv = np.zeros_like(i_local, dtype=np.int_)

    for count in prange(len(indices)):  # iterate over each given mesh element
        idx = indices[count]

        mass_elem = tetrahedral_volume_ref / 20 * det_b[idx] * value[count] * (np.ones((4, 4)) + np.eye(4))

        idx_node = mesh_elem2node[idx, :]
        idx_ijv = np.arange(16 * count, 16 * count + 16)
        i_local[idx_ijv] = np.column_stack((idx_node, idx_node, idx_node, idx_node)).flatten()
        j_local[idx_ijv] = np.hstack((idx_node, idx_node, idx_node, idx_node))
        v_local[idx_ijv] = mass_elem.flatten()
        indices_ijv[idx_ijv] = np.arange(16 * idx, 16 * idx + 16)

    i[indices_ijv] = i_local
    j[indices_ijv] = j_local
    v[indices_ijv] = v_local


def calc_mass_function_scalar(indices: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                              value: np.ndarray, det_b: np.ndarray, mesh_elem2node: np.ndarray,
                              weights: np.ndarray, local_coordinates: np.ndarray) -> None:
    """
    Calculate the elements of the curlcurl matrix if the material (reluctivity) is a function.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points per element
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    indices : np.ndarray
        (E,) array. Indices of mesh elements.
    i : np.ndarray
        (E,) array. Row-Index-Vector for sparse creation.
    j : np.ndarray
        (E,) array. Column-Index-Vector for sparse creation.
    v : np.ndarray
        (E,) array. Value-Vector for sparse creation.
    value: np.ndarray
        (E, N). Material value.
    det_b : np.ndarray
        (E,) array. Determinant of B-matrix.
    mesh_elem2node : np.ndarray
        (E,3) array. Indices of the nodes for each mesh element.
    weights : np.ndarray
        (N,) array. Weights for numerical integration.
    local_coordinates : np.ndarray
        (N,3) array. Local coordinates used for numerical integration.

    Returns
    -------
    None
    """
    i_local = np.zeros((16 * len(indices)), dtype=np.int_)
    j_local = np.zeros_like(i_local, dtype=np.int_)
    v_local = np.zeros_like(i_local, dtype=np.float_)
    indices_ijv = np.zeros_like(i_local, dtype=np.int_)

    ni = np.array([TetCartesianNodalShapeFunction.base1(local_coordinates[:, 0], local_coordinates[:, 1],
                                                        local_coordinates[:, 2]),
                   TetCartesianNodalShapeFunction.base2(local_coordinates[:, 0], local_coordinates[:, 1],
                                                        local_coordinates[:, 2]),
                   TetCartesianNodalShapeFunction.base3(local_coordinates[:, 0], local_coordinates[:, 1],
                                                        local_coordinates[:, 2]),
                   TetCartesianNodalShapeFunction.base4(local_coordinates[:, 0], local_coordinates[:, 1],
                                                        local_coordinates[:, 2])])

    for count in prange(len(indices)):  # iterate over each given mesh element
        idx = indices[count]

        mass_elem = det_b[idx] * ((value[count] * weights * ni) @ ni.T)

        idx_node = mesh_elem2node[idx, :]
        idx_ijv = np.arange(16 * count, 16 * count + 16)
        i_local[idx_ijv] = np.column_stack((idx_node, idx_node, idx_node, idx_node)).flatten()
        j_local[idx_ijv] = np.hstack((idx_node, idx_node, idx_node, idx_node))
        v_local[idx_ijv] = mass_elem.flatten()
        indices_ijv[idx_ijv] = np.arange(16 * idx, 16 * idx + 16)

    i[indices_ijv] = i_local
    j[indices_ijv] = j_local
    v[indices_ijv] = v_local


class TetCartesianNodalShapeFunction(NodalShapeFunction):
    """Class describing tetrahedral nodal shape functions on 3D cartesian meshs.

    Extends and implements the abstract class NodalShapeFunction to the case
    of first order tetrahedral nodal shape funcions, thus on a 3D cartesian mesh.
    Imagine classical 'hat functions' in 3D.

    Methods
    -------
    divgrad_operator:
        Compute the discrete version of the div-grad operator as a matrix.
    mass_matrix:
        Compute the discrete version of a mass term as a matrix.
    load_vector:
        Compute the discrete version of a node-wise excitation as a vector.
    gradient:
        Compute the gradient of the given nodal scalar field.
    robin_terms:

    See Also
    --------
    pyrit.shapefunction.ShapeFunction : Abstract super class
    pyrit.mesh : Abstract underlying Mesh object on whose entities the SFT
                  parameters are computed.
    """

    # highest_tetrahedral_integration_order = 3
    tetrahedral_volume_ref = 1 / 6  # Volume of unit tetrahedron
    dphi1 = np.array([[-1], [-1], [-1]])
    dphi2 = np.array([[1], [0], [0]])
    dphi3 = np.array([[0], [1], [0]])
    dphi4 = np.array([[0], [0], [1]])

    @staticmethod
    def issingleint(x) -> bool:
        """Workaround: check whether it is either an array with single integer or integer only."""
        if isinstance(x, (int, np.int_)):
            return True
        if isinstance(x, np.ndarray):
            if len(x.shape) == 1 and np.round(x) == x:
                return True
        return False

    @staticmethod
    def base1(x: float, y: float, z: float) -> float:
        """
        First basis function on reference element.

        Parameters
        ----------
        x : float
            x coordinate
        y : float
            y coordinate
        z : float
            z coordinate

        Returns
        -------
        float : value of function at given coordinates

        Notes
        -----
        Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf (German)
        """
        coord = np.array([x, y, z])
        if np.all(coord >= 0) and np.all(coord <= 1):
            return 1 - x - y - z
        return 0

    @staticmethod
    def base2(x: float, y: float, z: float) -> float:
        """
        Second basis function on reference element.

        Parameters
        ----------
        x : float
            x coordinate
        y : float
            y coordinate
        z : float
            z coordinate

        Returns
        -------
        float : value of function at given coordinates

        Notes
        -----
        Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf (German)
        """
        coord = np.array([x, y, z])
        if np.all(coord >= 0) and np.all(coord <= 1):
            return x
        return 0

    @staticmethod
    def base3(x: float, y: float, z: float) -> float:
        """
        Third basis function on reference element.

        Parameters
        ----------
        x : float
            x coordinate
        y : float
            y coordinate
        z : float
            z coordinate

        Returns
        -------
        float : value of function at given coordinates

        Notes
        -----
        Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf (German)
        """
        coord = np.array([x, y, z])
        if np.all(coord >= 0) and np.all(coord <= 1):
            return y
        return 0

    @staticmethod
    def base4(x: float, y: float, z: float) -> float:
        """
        Fourth basis function on reference element.

        Parameters
        ----------
        x : float
            x coordinate
        y : float
            y coordinate
        z : float
            z coordinate

        Returns
        -------
        float : value of function at given coordinates

        Notes
        -----
        Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf (German)
        """
        coord = np.array([x, y, z])
        if np.all(coord >= 0) and np.all(coord <= 1):
            return z
        return 0

    def __init__(self, mesh: TetMesh):
        super().__init__(mesh, dim=3, allocation_size=16 * mesh.num_elem)
        self.__determinant_b = np.zeros(mesh.num_elem)
        self.__inverse_b_transpose = np.zeros((3, 3, mesh.num_elem))
        self.__calc_transformation_matrices()

    def __calc_transformation_matrices(self):
        # get coordinates of nodes of tetrahedron
        p1 = self.mesh.node[self._mesh.elem2node[:, 0]]
        p2 = self.mesh.node[self._mesh.elem2node[:, 1]]
        p3 = self.mesh.node[self._mesh.elem2node[:, 2]]
        p4 = self.mesh.node[self._mesh.elem2node[:, 3]]
        # introduce shortcuts for difference vectors
        self.a = p2 - p1
        self.b = p3 - p1
        self.c = p4 - p1

        # for each element, define the local transformation matrix and calculate inverse and determinant
        for k in range(self._mesh.num_elem):
            matrix_b = np.c_[self.a[k, :], self.b[k, :], self.c[k, :]]
            self.__determinant_b[k] = np.abs(np.linalg.det(matrix_b))
            self.__inverse_b_transpose[:, :, k] = np.linalg.inv(matrix_b).T

    def _calc_matrix_constant_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: float, weights: np.ndarray, local_coordinates: np.ndarray,
                                     *args, evaluator: str = None, **kwargs):
        if matrix_type == "divgrad":
            calc_divgrad_constant_scalar(indices, i, j, v, value, self.__inverse_b_transpose, self.dphi1, self.dphi2,
                                         self.dphi3, self.dphi4, self.tetrahedral_volume_ref, self.__determinant_b,
                                         self.mesh.elem2node)
        elif matrix_type == "mass":
            calc_mass_constant_scalar(indices, i, j, v, value, self.tetrahedral_volume_ref, self.__determinant_b,
                                      self.mesh.elem2node)
        else:
            raise NotImplementedError(f"calc_matrix_constant_scalar is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_constant_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        raise NotImplementedError("Anisotrop materials are not implemented in TetCartesianNodalShapeFunction.")

    def _calc_matrix_scalar_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        if matrix_type == "divgrad":
            calc_divgrad_scalar_per_elem(indices, i, j, v, value, self.__inverse_b_transpose, self.dphi1, self.dphi2,
                                         self.dphi3, self.dphi4, self.tetrahedral_volume_ref, self.__determinant_b,
                                         self.mesh.elem2node)
        elif matrix_type == "mass":
            calc_mass_scalar_per_elem(indices, i, j, v, value, self.tetrahedral_volume_ref, self.__determinant_b,
                                      self.mesh.elem2node)
        else:
            raise NotImplementedError(f"calc_matrix_scalar_per_elem is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_tensor_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        raise NotImplementedError("Anisotrop materials are not implemented in TetCartesianNodalShapeFunction.")

    def _calc_matrix_function_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., float], weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        evaluation_points = calc_evaluation_points_element(local_coordinates, self.a, self.b, self.c,
                                                           self.mesh.node[self.mesh.elem2node[indices, 0]], indices)
        if evaluator == "eval_inhom_lin_iso":
            evaluations = eval_inhom_lin_iso(value, evaluation_points)
        else:
            raise NotImplementedError(f"Evaluator {evaluator} is not implemented.")

        if matrix_type == "divgrad":
            calc_divgrad_function_scalar(indices, i, j, v, evaluations, self.__inverse_b_transpose,
                                         self.dphi1, self.dphi2, self.dphi3, self.dphi4,
                                         self.__determinant_b, self.mesh.elem2node, weights)
        elif matrix_type == "mass":
            calc_mass_function_scalar(indices, i, j, v, evaluations, self.__determinant_b, self.mesh.elem2node,
                                      weights, local_coordinates)
        else:
            raise NotImplementedError(f"calc_matrix_function_scalar is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_function_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., np.ndarray], weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        raise NotImplementedError("Anisotrop materials are not implemented in TetCartesianNodalShapeFunction.")

    def divgrad_operator(self, *material: Union[Callable[..., float],
                                                np.ndarray, float, Tuple[Regions, Materials, MatProperty]],
                         integration_order: int = 1) -> coo_matrix:
        return self._matrix_routine("divgrad", *material, integration_order=integration_order)

    def mass_matrix(self, *material: Union[Callable[..., float],
                                           np.ndarray, float, Tuple['Regions', 'Materials', 'MatProperty']],
                    integration_order: int = 1) -> coo_matrix:
        return self._matrix_routine("mass", *material, integration_order=integration_order)

    def _calc_vec_elem_load_constant(self, load: Union[int, float], k: int) -> np.ndarray:
        """
        Calculates the local load vector for element k when the physical load / excitation is element-wise constant.

        Parameters
        ----------
        load : Union[int, float]
            Value of the load in element k
        k : int
            Index of element in mesh

        Returns
        -------
        vec_elem : np.ndarray
            (4,) array containing the load vector values for the four nodes of element k

        Notes
        -----
        For derivation, see info/Ansatzfunktionen/ansatzfunktionen
        """
        return self.tetrahedral_volume_ref / 4 * self.__determinant_b[k] * load

    def _calc_vec_elem_load_function(self, load: Callable[[float, float, float], float], k: int,
                                     weights: np.ndarray, int_points: np.ndarray) -> np.ndarray:
        """
        Calculates the local load vector for element k when the physical load / excitation is given as function.

        Parameters
        ----------
        load : Callable[[float, float, float], float]
            Function of physical load in element k in the format f(x,y,z)
        k : int
            Index of element in mesh
        weights : np.ndarray
            (m,) weights for numerical integration of order m over reference tetrahedron. See notes.
        int_points : np.ndarray
            (m, 3) integration points  of order m in reference tetrahedron. See notes.

        Returns
        -------
        vec_elem : np.ndarray
            (4,) array containing the load vector values for the four nodes of element k

        Notes
        -----
        For derivation, see info/Ansatzfunktionen/ansatzfunktionen.
        weights and int_points can be directly provided by pyrit/toolbox/quadraturetoolbox/gauss_tetrahedron
        """
        vec_elem = np.empty(4)
        matrix_b = np.c_[self.a[k, :], self.b[k, :], self.c[k, :]]  # local transformation matrix
        load_values = np.zeros(int_points.shape[0])
        for idx in range(int_points.shape[0]):  # iterate over integration points and find coordinates
            # in real mesh
            local_coord = int_points[idx, :].T
            real_coord = (matrix_b @ local_coord) + self.mesh.node[self.mesh.elem2node[k, 0]]
            load_values[idx] = load(real_coord[0], real_coord[1], real_coord[2])
        vec_elem[0] = self.__determinant_b[k] * np.sum(
            weights * load_values * self.base1(int_points[:, 0], int_points[:, 1], int_points[:, 2]))
        vec_elem[1] = self.__determinant_b[k] * np.sum(
            weights * load_values * self.base2(int_points[:, 0], int_points[:, 1], int_points[:, 2]))
        vec_elem[2] = self.__determinant_b[k] * np.sum(
            weights * load_values * self.base3(int_points[:, 0], int_points[:, 1], int_points[:, 2]))
        vec_elem[3] = self.__determinant_b[k] * np.sum(
            weights * load_values * self.base4(int_points[:, 0], int_points[:, 1], int_points[:, 2]))
        return vec_elem

    def load_vector(self, *load: Union[Callable[[float, float, float], float],
                                       np.ndarray, float, Tuple[Regions, Excitations]],
                    integration_order: int = 1) -> coo_matrix:
        # Input Check
        if len(load) == 1:
            case = self._process_arguments(*load, 2)
            load = load[0]
        elif len(load) == 2:
            if isinstance(load[-1], (int, np.int_)):
                raise ValueError("Wrong input format. Most likely you did not add 'integration-order=' in function "
                                 "call to distinguish it from load")
            case = self._process_arguments(load, 2)
        else:
            raise ValueError(f"Wrong number of inputs. Got {str(len(load))} but has to be 1 or 2.")
        if (integration_order != 1) and (case != "function"):
            warn("Integration order only makes a difference if a function is provided as of now.")
            # TODO: change this line once space-dependent excitations are implemented!
        weights, int_points = QuadratureToolbox.gauss_tetrahedron(int(integration_order))

        # setup of preliminary data vector
        v = np.zeros(self.mesh.num_node)
        # case distinction
        if case == 'array':
            for k in range(self.mesh.num_elem):
                node_ind = self.mesh.elem2node[k]
                v[node_ind] += self._calc_vec_elem_load_constant(load[k], k)
        elif case == 'number':
            for k in range(self.mesh.num_elem):
                node_ind = self.mesh.elem2node[k]
                v[node_ind] += self._calc_vec_elem_load_constant(load, k)
        elif case == 'function':
            for k in range(self.mesh.num_elem):
                node_ind = self.mesh.elem2node[k]
                vec_elem = self._calc_vec_elem_load_function(load, k, weights, int_points)
                # v[node_ind] += vec_elem
                v[node_ind[0]] += vec_elem[0]
                v[node_ind[1]] += vec_elem[1]
                v[node_ind[2]] += vec_elem[2]
                v[node_ind[3]] += vec_elem[3]
        elif case == 'tuple':
            regions, excitations = load
            for key in regions.get_keys():
                regi = regions.get_regi(key)
                exci_id = regi.exci
                if exci_id is None:  # if no excitation is assigned to region
                    continue
                exci = excitations.get_exci(exci_id)

                if regi.dim == 0:
                    indices = np.where(self.mesh.node2regi == regi.ID)[0]
                elif regi.dim == 1:
                    indices = np.where(self.mesh.edge2regi == regi.ID)[0]
                elif regi.dim == 2:
                    indices = np.where(self.mesh.face2regi == regi.ID)[0]
                elif regi.dim == 3:
                    indices = np.where(self.mesh.elem2regi == regi.ID)[0]
                else:
                    raise Exception(f"Internal error: unexpected region dimension {str(regi.dim)}. Should be one of "
                                    f"(0, 1, 2, 3). ")
                for k in indices:
                    if regi.dim == 0:
                        if exci.is_constant:
                            v[k] += exci.value  # note the += is necessary in case you have different dimensional
                            # excitations which affect the same node
                        else:  # currently only implemented for simple constant, homogeneous excitations!
                            raise NotImplementedError
                    elif regi.dim == 1:
                        if exci.is_constant:
                            node_ind = self.mesh.edge2node[k]
                            v[node_ind] += exci.value * (0.5 * self.mesh.edge_length[k])  # half the excitation of
                            # edge excitation for each node
                        else:
                            raise NotImplementedError
                    elif regi.dim == 2:
                        if exci.is_constant:
                            node_ind = self.mesh.face2node[k]
                            v[node_ind] += 1 / 3 * exci.value * self.mesh.face_area[k]
                        else:
                            raise NotImplementedError
                    elif regi.dim == 3:
                        if exci.is_constant:
                            node_ind = self.mesh.elem2node[k]
                            v[node_ind] += self._calc_vec_elem_load_constant(exci.value, k)
                        else:
                            raise NotImplementedError

        else:
            raise Exception(f"Severe Internal Error: returned case {case} does not match any of 'array', 'number', "
                            f"'function' or 'tuple'")

        lines = np.nonzero(v)[0]
        columns = np.zeros(len(lines), dtype=int)
        load = coo_matrix((v[lines], (lines, columns)), shape=(self.mesh.num_node, 1))

        return load

    def gradient(self, val2node: np.ndarray) -> np.ndarray:
        if not isinstance(val2node, np.ndarray):
            raise ValueError("Gradient is calculated only for array inputs.")
        if val2node.size != self.mesh.num_node:
            raise ValueError("Provided scalar field needs to provide one value for each node.")
        out = np.zeros((self.mesh.num_elem, 3))

        # Perform the simple gradient calculation element-wise
        for k in range(self.mesh.num_elem):
            invBTt = self.__inverse_b_transpose[:, :, k]
            # get node indices for this element
            n1, n2, n3, n4 = self.mesh.elem2node[k, :]
            out[k, :] = (val2node[n1] * invBTt @ self.dphi1
                         + val2node[n2] * invBTt @ self.dphi2
                         + val2node[n3] * invBTt @ self.dphi3
                         + val2node[n4] * invBTt @ self.dphi4).T

        return out

    def _calc_neumann_term_on_face_function(self, face_number, value, weights, int_points) -> np.ndarray:
        """
        Calculates the local neumann term vector for element k when a function is given.

        Parameters
        ----------
        face_number : int
            Index of element in mesh
        value : Callable[[float, float, float], float]
            Function of neumann term in element k in the format f(x,y,z)
        weights : np.ndarray
            (m,) weights for numerical integration of order m over reference triangle. See notes.
        int_points : np.ndarray
            (m, 2) integration points  of order m in reference triangle. See notes.

        Returns
        -------
        vec_elem : np.ndarray
            (3,) array containing the neumann vector values for the three nodes of face face_number

        Notes
        -----
        For derivation, see info/Ansatzfunktionen/ansatzfunktionen.
        weights and int_points can be directly provided by pyrit.toolbox.quadraturetoolbox.gauss_triangle
        """
        vec_elem = np.zeros(4)
        face_nodes = self.mesh.face2node[face_number, :]
        matrix_b = np.c_[self.mesh.node[face_nodes[1], :] - self.mesh.node[face_nodes[0], :],
                         self.mesh.node[face_nodes[2], :] - self.mesh.node[face_nodes[0], :],
                         np.ones((3,))]  # local transformation matrix
        nm_values = np.zeros(weights.size)
        for idx in range(int_points.shape[0]):
            # iterate over integration points and find coordinates in real mesh
            local_coord = np.r_[int_points[idx, :], 0]
            real_coord = (matrix_b @ local_coord) + self.mesh.node[face_nodes[0]]
            nm_values[idx] = value(*real_coord)

        dummyz = np.zeros(int_points[:, 0].shape)
        vec_elem[0] = np.sum(weights * nm_values * self.base1(int_points[:, 0], int_points[:, 1], dummyz))
        vec_elem[1] = np.sum(weights * nm_values * self.base2(int_points[:, 0], int_points[:, 1], dummyz))
        vec_elem[2] = np.sum(weights * nm_values * self.base3(int_points[:, 0], int_points[:, 1], dummyz))
        vec_elem = vec_elem * np.abs(np.linalg.det(matrix_b))
        return vec_elem

    def _calc_neumann_term_on_face_constant(self, face_number, value) -> np.ndarray:
        return np.ones(3) * self.mesh.face_area[face_number] * value / 3

    # pylint: disable=cell-var-from-loop
    def neumann_term(self, *args: Union[Tuple[Regions, BdryCond],
                                        Tuple[np.ndarray, Union[Callable[..., float], np.ndarray, float]]],
                     integration_order: int = 1) -> coo_matrix:
        # note complete similarity to implementation in TriCartesianNodalShapeFunction
        flag_regions, flag_value = self._process_neumann(*args)

        weights, int_points = QuadratureToolbox.gauss_triangle(int(integration_order))

        neumann_vector = np.zeros((self.mesh.num_node, 1))
        if flag_regions:  # We have regions and boundary_conditions
            regions, boundary_conditions = args
            regions: Regions
            boundary_conditions: BdryCond
            for key in boundary_conditions.get_ids():
                bc = boundary_conditions.get_bc(key)
                if isinstance(bc, BCNeumann):
                    bc_value = bc.value
                    indices_regions = regions.find_regions_of_boundary_condition(bc.ID)
                    indices_faces = np.empty(0)  # Indices of all faces that are on the boundary bc
                    for region_id in indices_regions:
                        indices_faces = np.r_[indices_faces, np.where(self.mesh.face2regi == region_id)[0]]
                    indices_faces = indices_faces.astype('int')
                    if not bc.is_homogeneous:
                        for face in indices_faces:
                            face_values = self._calc_neumann_term_on_face_function(face, bc_value, weights, int_points)
                            neumann_vector[self.mesh.face2node[face, 0], 0] += face_values[0]
                            neumann_vector[self.mesh.face2node[face, 1], 0] += face_values[1]
                            neumann_vector[self.mesh.face2node[face, 2], 0] += face_values[2]
                    if bc.is_constant:
                        for k, face in enumerate(indices_faces):
                            if isinstance(bc_value, np.ndarray):
                                face_values = self._calc_neumann_term_on_face_constant(face, bc_value[k])
                            elif isinstance(bc_value, (float, int, complex)):
                                face_values = self._calc_neumann_term_on_face_constant(face, bc_value)
                            else:
                                raise Exception("Type not supported")
                            neumann_vector[self.mesh.face2node[face, 0], 0] += face_values[0]
                            neumann_vector[self.mesh.face2node[face, 1], 0] += face_values[1]
                            neumann_vector[self.mesh.face2node[face, 2], 0] += face_values[2]
        else:  # We have indices and value
            indices, value = args
            for k, index in enumerate(indices):
                if flag_value == 'callable':
                    node_test = self.mesh.node[self.mesh.face2node[index, 0], :]
                    if isinstance(value(*node_test), np.ndarray):
                        # noinspection PyCallingNonCallable
                        face_values = self._calc_neumann_term_on_face_function(index, lambda x, y, z: value(x, y, z)[k],
                                                                               weights, int_points)
                    else:
                        face_values = self._calc_neumann_term_on_face_function(index, value, weights, int_points)
                elif flag_value == 'array':
                    face_values = self._calc_neumann_term_on_face_constant(index, value[k])
                elif flag_value == 'value':
                    face_values = self._calc_neumann_term_on_face_constant(index, value)
                # noinspection PyUnboundLocalVariable
                neumann_vector[self.mesh.face2node[index, 0], 0] += face_values[0]
                neumann_vector[self.mesh.face2node[index, 1], 0] += face_values[1]
                neumann_vector[self.mesh.face2node[index, 2], 0] += face_values[2]
        return coo_matrix(neumann_vector)

    def _calc_robin_term_function(self, face_number, value, weights, int_points):
        """Calculates the local Robin matrix vector for a surface when a function is given.

        Parameters
        ----------
        face_number : int
            Index of element in mesh
        value : Callable[[float, float, float], float]
            Function of Robin alpha/beta term in element k in the format f(x,y,z)
        weights : np.ndarray
            (m,) weights for numerical integration of order m over reference triangle. See notes.
        int_points : np.ndarray
            (m, 2) integration points  of order m in reference triangle. See notes.


        Returns
        -------
        vec_elem : np.ndarray
            (9,) array containing the Robin matrix values for the three nodes of face face_number

        Notes
        -----
        Compare the 2D TriCartesian Mass matrix.
        weights and int_points can be directly provided by pyrit.toolbox.quadraturetoolbox.gauss_triangle
        """
        face_nodes = self.mesh.face2node[face_number, :]
        matrix_b = np.c_[self.mesh.node[face_nodes[1], :] - self.mesh.node[face_nodes[0], :],
                         self.mesh.node[face_nodes[2], :] - self.mesh.node[face_nodes[0], :],
                         np.ones((3,))]  # local transformation matrix
        # compare TriCart Mass-Matrix
        integral = np.zeros((3, 3))
        for kk, weight in enumerate(weights):
            m_values = np.array(
                [1 - int_points[kk, 0] - int_points[kk, 1], int_points[kk, 0],
                 int_points[kk, 1]])
            local_coord = np.array(int_points[kk, :], ndmin=2).T
            tmp = (matrix_b[:, 0:2] @ local_coord).flatten() + self.mesh.node[face_nodes[0]]
            # noinspection PyCallingNonCallable
            evaluation = value(x=tmp[0], y=tmp[1], z=tmp[2])
            integral = integral + weight * evaluation * np.outer(m_values, m_values)
        return np.abs(np.linalg.det(matrix_b)) * integral.reshape(9)

    def _calc_robin_term_constant(self, face_number, value):
        """Calculates the local Robin matrix vector for a surface when a constant is given.

        Parameters
        ----------
        face_number : int
            Index of element in mesh
        value : float
            Constant value alpha / beta

        Returns
        -------
        vec_elem : np.ndarray
            (9,) array containing the Robin matrix values for the three nodes of face face_number

        Notes
        -----
        Compare the 2D TriCartesian Mass matrix.
        """
        return self.mesh.face_area[face_number] * value / 12 * np.array([2, 1, 1, 1, 2, 1, 1, 1, 2])

    # pylint: disable=unnecessary-lambda-assignment
    def robin_terms(self, regions: Regions, boundary_condition: 'BdryCond', integration_order: int = 1) -> \
            Tuple[coo_matrix, coo_matrix]:

        # simple input check: arguments are of correct class
        if not (isinstance(regions, Regions) and isinstance(boundary_condition, BdryCond)):
            raise ValueError("Input arguments have wrong types. Did you mix up the order?")

        bc_keys = boundary_condition.get_ids()

        values_list = []
        lines_list = []
        columns_list = []
        robin_vector = coo_matrix((self.mesh.num_node, 1))
        weights, int_points = QuadratureToolbox.gauss_triangle(int(integration_order))

        def check_valid_function(arg, numargs):
            if callable(arg):
                if len(signature(arg).parameters) != numargs:
                    raise ValueError(f"The function has not the expected number of {str(numargs)} "
                                     f"arguments")

        for bc_key in bc_keys:
            bc = boundary_condition.get_bc(bc_key)
            if not isinstance(bc, BCRobin):
                continue

            # Calculate the matrix
            indices_regions = regions.find_regions_of_boundary_condition(bc.ID)
            indices_faces = np.empty(0)  # Indices of all faces that are on the boundary bc
            for region_id in indices_regions:
                indices_faces = np.r_[indices_faces, np.where(self.mesh.face2regi == region_id)[0]]
            indices_faces = indices_faces.astype('int')
            robin_values_per_bc = np.zeros(9 * len(indices_faces))
            robin_lines_per_bc = np.zeros(9 * len(indices_faces), dtype=int)
            robin_columns_per_bc = np.zeros(9 * len(indices_faces), dtype=int)

            coef_dir = bc.coefficient_dirichlet
            coef_neum = bc.coefficient_neumann
            if np.any(coef_neum == 0):
                raise ValueError("Neumann coefficient must not be zero.")
            value_bc = bc.value
            if callable(coef_dir):
                check_valid_function(coef_dir, 3)
                if callable(coef_neum):
                    check_valid_function(coef_neum, 3)
                    # Both coefficients are functions
                    for k, face in enumerate(indices_faces):
                        idx_node = np.kron(np.ones((3, 1), dtype=int), self.mesh.face2node[face])
                        tmp_values = self._calc_robin_term_function(face,
                                                                    lambda x, y, z: coef_dir(x, y, z) / coef_neum(x, y,
                                                                                                                  z),
                                                                    weights, int_points)
                        robin_values_per_bc[9 * k:9 * (k + 1)] = tmp_values
                        robin_lines_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node, 9)
                        robin_columns_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node.T, 9)
                else:
                    # Only dirichlet coefficient is a function
                    if isinstance(coef_neum, np.ndarray):
                        if len(indices_faces) != coef_neum.size:
                            raise ValueError("Provided array-valued coefficient does not have the same number of "
                                             "elements as there are surface elements for the Robin BC.")
                    for k, face in enumerate(indices_faces):
                        idx_node = np.kron(np.ones((3, 1), dtype=int), self.mesh.face2node[face])
                        if isinstance(coef_neum, np.ndarray):
                            tmp_values = self._calc_robin_term_function(face,
                                                                        lambda x, y, z: coef_dir(x, y, z) / coef_neum[
                                                                            k], weights, int_points)
                        elif isinstance(coef_neum, (float, int)):
                            tmp_values = self._calc_robin_term_function(face,
                                                                        lambda x, y, z: coef_dir(x, y, z) / coef_neum,
                                                                        weights, int_points)
                        else:
                            raise Exception("Type of coef_neum is not supported")
                        robin_values_per_bc[9 * k:9 * (k + 1)] = tmp_values
                        robin_lines_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node, 9)
                        robin_columns_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node.T, 9)
            else:
                if callable(coef_neum):
                    check_valid_function(coef_neum, 3)
                    if isinstance(coef_dir, np.ndarray):
                        if len(indices_faces) != coef_dir.size:
                            raise ValueError("Provided array-valued coefficient does not have the same number of "
                                             "elements as there are surface elements for the Robin BC.")
                    # Only coef_neum is a function
                    for k, face in enumerate(indices_faces):
                        idx_node = np.kron(np.ones((3, 1), dtype=int), self.mesh.face2node[face])
                        if isinstance(coef_dir, np.ndarray):
                            tmp_values = self._calc_robin_term_function(face,
                                                                        lambda x, y, z: coef_dir[k] / coef_neum(x, y,
                                                                                                                z),
                                                                        weights, int_points)
                        elif isinstance(coef_dir, (float, int)):
                            tmp_values = self._calc_robin_term_function(face,
                                                                        lambda x, y, z: coef_dir / coef_neum(x, y, z),
                                                                        weights, int_points)
                        else:
                            raise Exception("Type of coef_dir is not supported")
                        robin_values_per_bc[9 * k:9 * (k + 1)] = tmp_values
                        robin_lines_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node, 9)
                        robin_columns_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node.T, 9)
                else:
                    # Both are no functions
                    value = coef_dir / coef_neum
                    if isinstance(value, np.ndarray):
                        if len(indices_faces) != value.size:
                            raise ValueError("Provided array-valued coefficient does not have the same number of "
                                             "elements as there are surface elements for the Robin BC.")
                    for k, face in enumerate(indices_faces):
                        idx_node = np.kron(np.ones((3, 1), dtype=int), self.mesh.face2node[face])
                        if isinstance(value, np.ndarray):
                            tmp_values = self._calc_robin_term_constant(face, value[k])
                        elif isinstance(value, (float, int)):
                            tmp_values = self._calc_robin_term_constant(face, value)
                        else:
                            raise Exception("Type of value not supported")
                        robin_values_per_bc[9 * k:9 * (k + 1)] = tmp_values
                        robin_lines_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node, 9)
                        robin_columns_per_bc[9 * k:9 * (k + 1)] = np.reshape(idx_node.T, 9)
            values_list.append(robin_values_per_bc)
            lines_list.append(robin_lines_per_bc)
            columns_list.append(robin_columns_per_bc)

            # Calculate the vector
            if not bc.is_homogeneous:
                if callable(coef_neum):
                    check_valid_function(coef_neum, 3)
                    check_valid_function(value_bc, 3)
                    # Both are functions
                    if callable(value_bc):
                        fun = lambda x, y, z: value_bc(x, y, z) / coef_neum(x, y, z)
                    else:
                        fun = lambda x, y, z: value_bc / coef_neum(x, y, z)
                    robin_vector = robin_vector + self.neumann_term(indices_faces, fun,
                                                                    integration_order=integration_order)
                else:
                    # only value_bc is a function
                    check_valid_function(value_bc, 3)
                    if callable(value_bc):
                        fun = lambda x, y, z: value_bc(x, y, z) / coef_neum
                    else:
                        fun = value_bc / coef_neum
                    robin_vector = robin_vector + self.neumann_term(indices_faces, fun,
                                                                    integration_order=integration_order)
            else:
                if callable(coef_neum):
                    check_valid_function(coef_neum, 3)
                    # Only Neumann coefficient is a function
                    if isinstance(value_bc, np.ndarray):
                        robin_vector = robin_vector + self.neumann_term(indices_faces,
                                                                        lambda x, y, z: value_bc / coef_neum(x, y, z),
                                                                        integration_order=integration_order)
                    elif isinstance(value_bc, (float, int)):
                        robin_vector = robin_vector + self.neumann_term(indices_faces,
                                                                        lambda x, y, z: value_bc / coef_neum(x, y, z),
                                                                        integration_order=integration_order)
                    else:
                        raise Exception("Type of value not supported")

                else:
                    value_tmp = value_bc / coef_neum
                    if isinstance(value_tmp, np.ndarray):
                        robin_vector = robin_vector + self.neumann_term(indices_faces, value_tmp,
                                                                        integration_order=integration_order)
                    elif isinstance(value_tmp, (float, int)):
                        robin_vector = robin_vector + self.neumann_term(indices_faces, value_tmp,
                                                                        integration_order=integration_order)
                    else:
                        raise Exception("Type not supported.")

        if len(values_list) != 0:
            robin_matrix = coo_matrix(
                (np.concatenate(values_list), (np.concatenate(lines_list), np.concatenate(columns_list))),
                shape=(self.mesh.num_node, self.mesh.num_node))
        else:
            robin_matrix = coo_matrix((self.mesh.num_node, self.mesh.num_node)).tocoo()
        return robin_matrix, robin_vector.tocoo()
