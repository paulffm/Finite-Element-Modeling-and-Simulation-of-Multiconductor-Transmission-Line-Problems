# coding=utf-8
"""
File containing the class TriAxisymmetricNodalShapeFunction.

Numba is used with a lot functions to get a good performance. The data types of the arguments of the numa functions
can be given in the decorator. If one of these functions is called with another data type it does not work anymore.
Because different data types of some arguments should be supported the data types are not indicated in the decorator.
Anyway, the code is written such that the data types of the arguments are as equal as possible for different user input
parameters. So the number of compilations with numba should be at a minimum.

.. sectionauthor:: bundschuh, ruppert
"""
# pylint: disable=duplicate-code

from typing import Union, Callable, Tuple, NoReturn, TYPE_CHECKING

import numpy as np
from numpy import ndarray
from scipy.sparse import coo_matrix, csr_matrix
from numba import prange

from pyrit.region import Regions
from pyrit.bdrycond import BdryCond, BCNeumann, BCRobin

from pyrit.toolbox.QuadratureToolbox import gauss_triangle, gauss
from pyrit.toolbox.NumbaToolbox import njit_cache

from . import NodalShapeFunction

if TYPE_CHECKING:
    from pyrit.mesh import AxiMesh
    from pyrit.material import MatProperty, Materials
    from pyrit.excitation import Excitations


# region General functions
@njit_cache(parallel=True)
def calc_evaluation_points_element(local_coordinates: np.ndarray, transform_coefficients: np.ndarray,
                                   elements: np.ndarray) -> np.ndarray:
    """Calculates the coordinates (rho,z) needed for the evaluation of the numerical integration on an element.

    Returns a three dimensional array with the evaluation coordinates for all elements specified. The return values is
    a (E,N,2) array. It is `out[e,n,0]` and `out[e,n,1]` the rho and z coordinate on element e and evaluation point
    n, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        ======  =======

    Parameters
    ----------
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    transform_coefficients : np.ndarray
        (T,3,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elements : np.ndarray
        (E,) array. Array with the indices of elements the evalutaion points should calculated for.

    Returns
    -------
    out : np.ndarray
        (E,N,2) array.
    """
    num_elements = len(elements)
    num_coordinates = local_coordinates.shape[0]
    coords = np.empty((num_elements, num_coordinates, 2))
    ones = np.ones(num_coordinates)
    local_cds = np.column_stack((local_coordinates, ones))
    for k in prange(num_elements):
        element = elements[k]
        tc = transform_coefficients[element]  # (2,3) array
        tct = tc.T
        coords[k] = local_cds @ tct
    return coords


@njit_cache(parallel=True)
def calc_evaluation_points_edge(local_coordinates: np.ndarray, nodes: np.ndarray,
                                edges: np.ndarray, edge2node: np.ndarray) -> np.ndarray:
    """Calculates the coordinates (s,z) needed for the evaluation of the numerical integration on an edge.

    Returns a three dimensional array with the evaluation coordinates for all edges specified. The return values is
    a (E,N,2) array. It is `out[e,n,0]` and `out[e,n,1]` the s and z coordinate on edge e and evaluation point n,
    respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of edges in `edges`
        T       Number of elements in the mesh
        ======  =======

    Parameters
    ----------
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    node_transformed : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edges : np.ndarray
        (E,) array. Array with the indices of edges the evalutaion points should calculated for.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    out : np.ndarray
        (E,N,2) array.
    """
    num_edges = len(edges)
    num_coordinates = len(local_coordinates)
    coords = np.empty((num_edges, num_coordinates, 2))

    for k in prange(num_edges):
        edge = edges[k]
        for kk in range(num_coordinates):
            coords[k, kk] = nodes[edge2node[edge, 0]] + local_coordinates[kk] * \
                (nodes[edge2node[edge, 1]] - nodes[edge2node[edge, 0]])

    return coords


def eval_hom_nonlin_iso(fun: Callable[[np.ndarray], np.ndarray], elements: np.ndarray):
    """Evaluates the function `fun` on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[int], float]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

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
    # noinspection PyArgumentList
    return fun(element=elements)


def eval_hom_nonlin_aniso(fun: Callable[[np.ndarray], np.ndarray], elements: np.ndarray):
    """Evaluates the function `fun` on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[int], np.ndarray]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

    Returns
    -------
    evaluations : np.ndarray
        (E,N,2,2) array. For every element the N evaluations of the material function at the evaluation points for the
        numerical integration.
    """
    # noinspection PyArgumentList
    return fun(element=elements)


def eval_inhom_lin_iso(fun: Callable[[float, float], float], evaluation_points: np.ndarray):
    """Evaluates the function `fun` for numerical integration on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float], float]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

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
    evaluations = np.vectorize(fun)(evaluation_points[:, :, 0], evaluation_points[:, :, 1])
    # evaluations.astype(float, copy=False)  #change: why make it float?
    return evaluations


def eval_inhom_lin_aniso(fun: Callable[[float, float], np.ndarray], evaluation_points: np.ndarray):
    """Evaluates the function `fun` for numerical integration on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float], np.ndarray]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

    Returns
    -------
    evaluations : np.ndarray
        (E,N,2,2) array. For every element the N evaluations of the material function at the evaluation points for the
        numerical integration.
    """
    vfun = np.vectorize(fun, signature='(),()->(i,i)')
    evaluations = vfun(evaluation_points[:, :, 0], evaluation_points[:, :, 1])
    # evaluations.astype(float, copy=False)  #change: why make it float?
    return evaluations


def eval_inhom_nonlin_iso(fun: Callable[[float, float, int], float], evaluation_points: np.ndarray,
                          elements):
    """Evaluates the function `fun` for numerical integration on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float], float]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.
    elements : np.ndarray
        Array of the considered elements.

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
    indices = np.kron(np.ones((evaluation_points.shape[1], 1), dtype=int), elements).transpose()
    evaluations = np.vectorize(fun)(evaluation_points[:, :, 0], evaluation_points[:, :, 1], indices)
    # evaluations.astype(float, copy=False)  #change: why make it float?
    return evaluations


# todo:docstring
def eval_inhom_nonlin_aniso(fun: Callable[[float, float, int], np.ndarray], evaluation_points: np.ndarray,
                            elements):
    """Evaluates the function `fun` for numerical integration on elements.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified elements.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of elements in `elements`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float, int], np.ndarray]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.

    Returns
    -------
    evaluations : np.ndarray
        (E,N,2,2) array. For every element the N evaluations of the material function at the evaluation points for the
        numerical integration.
    """
    indices = np.kron(np.ones((evaluation_points.shape[1], 1), dtype=int), elements).transpose()
    vfun = np.vectorize(fun, signature='(),(),()->(i,i)')
    evaluations = vfun(evaluation_points[:, :, 0], evaluation_points[:, :, 1], indices)
    # evaluations.astype(float, copy=False)  #change: why make it float?
    return evaluations


def evaluate_for_num_int_edge(fun: Callable[[float, float], float], evaluation_points: np.ndarray):
    """Evaluates the function `fun` for numerical integration on edges.

    Evaluates `fun` at the evaluation points for the numerical integration for all specified edges.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        N       Number of evaluation points
        E       Number of edges in `edges`
        ======  =======

    Parameters
    ----------
    fun : Callable[[float, float], float]
        The function to evaluate.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_edge`.

    Returns
    -------
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration.

    Notes
    -----
    This function has to be called before matrix elements, where a numerical integration of a arbitrary funtion is
    necessary, can be computed.
    """
    evaluations = np.vectorize(fun)(evaluation_points[:, :, 0], evaluation_points[:, :, 1])
    evaluations.astype(float, copy=False)
    return evaluations


# endregion


# region Functions for divgrad


@njit_cache(parallel=True)
def calc_divgrad_constant_scalar(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: Union[int, float], b: np.ndarray, c: np.ndarray,
                                 bmat_11: np.ndarray, bmat_12: np.ndarray,
                                 r1: np.ndarray, determinant_b: np.ndarray,
                                 elem2node) -> NoReturn:
    """Calculates elements for the divgrad_operator with a constant material.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : Union[int, float]
        The value of the material.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_11: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9), np.int_)
    j_local = np.zeros((num_elements, 9), np.int_)
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        b_mat = np.outer(b[element], b[element])
        c_mat = np.outer(c[element], c[element])

        mat_elem = 2 / determinant_b[element] * value * np.pi * \
                   ((bmat_11[element] + bmat_12[element]) / 6 + 0.5 * r1[element]) * (c_mat + b_mat)
        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


# todo:docstring
@njit_cache(parallel=True)
def calc_divgrad_constant_tensor(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: np.ndarray, b: np.ndarray, c: np.ndarray,
                                 bmat_11: np.ndarray, bmat_12: np.ndarray,
                                 r1: np.ndarray, determinant_b: np.ndarray,
                                 elem2node) -> NoReturn:
    """
    Calculates elements for the divgrad_operator with a constant material tensor.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : np.ndarray
        (2,2) The value of the material.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_11: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9), np.int_)
    j_local = np.zeros((num_elements, 9), np.int_)
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        b_mat = np.outer(b[element], b[element])
        c_mat = np.outer(c[element], c[element])
        bc_mat = np.outer(b[element], c[element])
        cb_mat = np.outer(c[element], b[element])

        mat_elem = 2 / determinant_b[element] * np.pi * \
                   ((bmat_11[element] + bmat_12[element]) / 6 + 0.5 * r1[element]) * \
                   (b_mat * value[0, 0] + cb_mat * value[1, 0] + bc_mat * value[0, 1] + c_mat * value[1, 1])
        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_divgrad_scalar_per_elem(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: np.ndarray, b: np.ndarray, c: np.ndarray,
                                 bmat_11: np.ndarray, bmat_12: np.ndarray,
                                 r1: np.ndarray, determinant_b: np.ndarray,
                                 elem2node):
    """Calculates elements for the divgrad_operator with one constant value of the material per element.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : np.ndarray
        (E,) array. One material value per element.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_11 : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12 : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1 : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9))
    j_local = np.zeros((num_elements, 9))
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        b_mat = np.outer(b[element], b[element])
        c_mat = np.outer(c[element], c[element])

        mat_elem = 2 / determinant_b[element] * value[k] * np.pi * \
                   ((bmat_11[element] + bmat_12[element]) / 6 + 0.5 * r1[element]) * \
                   (c_mat + b_mat)

        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


# todo:docstring
@njit_cache(parallel=True)
def calc_divgrad_tensor_per_elem(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 value: np.ndarray, b: np.ndarray, c: np.ndarray,
                                 bmat_11: np.ndarray, bmat_12: np.ndarray,
                                 r1: np.ndarray, determinant_b: np.ndarray,
                                 elem2node) -> NoReturn:
    """
    Calculates elements for the divgrad_operator with one constant material tensor per element.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : np.ndarray
        (T,2,2) The value of the material.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_11: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1: np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9), np.int_)
    j_local = np.zeros((num_elements, 9), np.int_)
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        b_mat = np.outer(b[element], b[element])
        c_mat = np.outer(c[element], c[element])
        bc_mat = np.outer(b[element], c[element])
        cb_mat = np.outer(c[element], b[element])

        mat_elem = 2 / determinant_b[element] * np.pi * \
                   ((bmat_11[element] + bmat_12[element]) / 6 + 0.5 * r1[element]) * \
                   (b_mat * value[k, 0, 0] + cb_mat * value[k, 1, 0] + bc_mat * value[k, 0, 1] + c_mat * value[k, 1, 1])
        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_divgrad_function_scalar(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 weights: np.ndarray, evaluation_points: np.ndarray,
                                 evaluations: np.ndarray, b: np.ndarray, c: np.ndarray,
                                 determinant_b: np.ndarray, elem2node: np.ndarray):
    """Calculates elements for the divgrad_operator with a function as material.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_element`.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9))
    j_local = np.zeros((num_elements, 9))
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        index = elements[k]

        b_mat = np.outer(b[index], b[index])
        c_mat = np.outer(c[index], c[index])

        integral = np.dot(weights, evaluations[k] * evaluation_points[k, :, 0])
        mat_elem = 2 / determinant_b[index] * np.pi * (c_mat + b_mat) * integral
        idx_node = elem2node[index, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * index, 9 * index + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)  # todo:docstring
def calc_divgrad_function_tensor(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                 weights: np.ndarray, evaluation_points: np.ndarray,
                                 evaluations: np.ndarray, b: np.ndarray, c: np.ndarray,
                                 determinant_b: np.ndarray, elem2node: np.ndarray):
    """Calculates elements for the divgrad_operator with a function as material.

    Calculates the contribution of the elements specified in `elements` to the divgrad matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_element`.
    b : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9))
    j_local = np.zeros((num_elements, 9))
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        b_mat = np.outer(b[element], b[element])
        c_mat = np.outer(c[element], c[element])
        bc_mat = np.outer(b[element], c[element])
        # cb_mat = np.outer(c[element], b[element])
        # ToDo: Check if cb_mat can always be replaced by bc_mat

        first_integral = np.dot(weights, evaluations[k, :, 0, 0] * evaluation_points[k, :, 0])
        second_integral = np.dot(weights, evaluations[k, :, 0, 1] * evaluation_points[k, :, 0])
        third_integral = np.dot(weights, evaluations[k, :, 1, 0] * evaluation_points[k, :, 0])
        fourth_integral = np.dot(weights, evaluations[k, :, 1, 1] * evaluation_points[k, :, 0])

        sum_integrals = b_mat * first_integral + bc_mat * second_integral - \
                        bc_mat * third_integral + c_mat * fourth_integral

        mat_elem = 2 / determinant_b[element] * np.pi * sum_integrals
        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


# endregion


# region Functions for mass
@njit_cache(parallel=True)
def calc_mass_constant_scalar(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                              value: Union[int, float],
                              bmat_11: np.ndarray, bmat_12: np.ndarray,
                              r1: np.ndarray, determinant_b: np.ndarray,
                              elem2node):
    """Calculates elements for the mass matrix with one constant value of the material per element.

    Calculates the contribution of the elements specified in `elements` to the mass matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : np.ndarray
        (E,) array. One value per element.
    bmat_11 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9))
    j_local = np.zeros((num_elements, 9))
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        mat_elem = determinant_b[element] * value * np.pi * \
                   (bmat_11[element] / 60 * np.array([[2, 2, 1], [2, 6, 2], [1, 2, 2]]) +
                    bmat_12[element] / 60 * np.array([[2, 1, 2], [1, 2, 2], [2, 2, 6]]) +
                    r1[element] / 12 * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]))

        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


def calc_mass_scalar_per_elem(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                              value: np.ndarray,
                              bmat_11: np.ndarray, bmat_12: np.ndarray,
                              r1: np.ndarray, determinant_b: np.ndarray,
                              elem2node):
    """Calculates elements for the mass matrix with one constant value of the material per element.

    Calculates the contribution of the elements specified in `elements` to the mass matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.divgrad_operator`.
    value : np.ndarray
        (E,) array. One value per element.
    bmat_11 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1 : np.ndarray
        (T,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9))
    j_local = np.zeros((num_elements, 9))
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        mat_elem = determinant_b[element] * value[k] * np.pi * \
                   (bmat_11[element] / 60 * np.array([[2, 2, 1], [2, 6, 2], [1, 2, 2]]) +
                    bmat_12[element] / 60 * np.array([[2, 1, 2], [1, 2, 2], [2, 2, 6]]) +
                    r1[element] / 12 * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]))

        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_mass_function_scalar(elements: np.ndarray, i: np.ndarray, j: np.ndarray, v: np.ndarray,
                                    weights: np.ndarray, local_coordinates: np.ndarray,
                                    evaluation_points: np.ndarray, evaluations: np.ndarray,
                                    determinant_b: np.ndarray, elem2node: np.ndarray):
    """Calculates elements for the mass matrix with a function as material.

    Calculates the contribution of the elements specified in `elements` to the mass matrix. The calculated values
    and line and column information are written to v, i and j, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements : np.ndarray
        (E,) array. The indices of elements that should be considered
    i : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.mass_matrix`.
    j : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.mass_matrix`.
    v : np.ndarray
        (K,) array. See the variable in :py:func:`TriAxisymmetricNodalShapeFunction.mass_matrix`.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_element`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_ijv = np.zeros((num_elements, 9), np.int_)
    i_local = np.zeros((num_elements, 9), np.int_)
    j_local = np.zeros((num_elements, 9), np.int_)
    v_local = np.zeros((num_elements, 9))
    for k in prange(num_elements):
        element = elements[k]

        integral = np.zeros((3, 3))
        for kk, weight in enumerate(weights):
            m_values = np.array([1 - local_coordinates[kk, 0] - local_coordinates[kk, 1], local_coordinates[kk, 0],
                                 local_coordinates[kk, 1]])
            integral = integral + weight * evaluations[k, kk] * evaluation_points[k, kk, 0] * np.outer(m_values,
                                                                                                       m_values)

        mat_elem = determinant_b[element] * 2 * np.pi * integral
        idx_node = elem2node[element, :]

        tmp = np.stack((idx_node, idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node, idx_node))

        indices_ijv[k] = np.arange(9 * element, 9 * element + 9, 1)
        i_local[k] = np.reshape(tmp, (9,))
        j_local[k] = np.reshape(tmpt, (9,))
        v_local[k] = np.reshape(mat_elem, (9,))

    # Write results to i, j, v
    idx = indices_ijv.flatten()
    i[idx] = i_local.flatten()
    j[idx] = j_local.flatten()
    v[idx] = v_local.flatten()


# endregion


# region Functions for load vector

@njit_cache(parallel=True)
def calc_load_elems_constant(elements: np.ndarray, i: np.ndarray, v: np.ndarray,
                             value: Union[int, float], bmat_11: np.ndarray, bmat_12: np.ndarray,
                             r1: np.ndarray, determinant_b: np.ndarray,
                             elem2node) -> NoReturn:
    """Calculates elements for the load vector with a constant load.

    Calculates the contribution of the elements specified in `elements` to the load vector. The calculated values and
    line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements :  np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable `i_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See the variable `v_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    value :  Union[int, float]
        The value of the material.
    bmat_11 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)
    indices_v = np.zeros((num_elements, 3), np.int_)
    i_local = np.zeros((num_elements, 3))
    v_local = np.zeros((num_elements, 3))
    for k in prange(num_elements):
        element = elements[k]
        vec_elem = value / 12 * np.pi * determinant_b[element] * (bmat_11[element] * np.array([1, 2, 1]) +
                                                                  bmat_12[element] * np.array([1, 1, 2]) +
                                                                  r1[element] * np.array([4, 4, 4]))
        idx_node = elem2node[element, :]

        indices_v[k] = np.arange(3 * element, 3 * element + 3, 1)
        i_local[k] = idx_node
        v_local[k] = vec_elem

    # Write results to v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_load_elems_vector(elements: np.ndarray, i: np.ndarray, v: np.ndarray,
                           value: Union[int, float], bmat_11: np.ndarray, bmat_12: np.ndarray,
                           r1: np.ndarray, determinant_b: np.ndarray,
                           elem2node) -> NoReturn:
    """Calculates elements for the load vector with a constant load.

    Calculates the contribution of the elements specified in `elements` to the load vector. The calculated values and
    line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements :  np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable `i_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See the variable `v_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    value :  Union[int, float]
        The value of the material.
    bmat_11 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    bmat_12 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    r1 : np.ndarray
        (N,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)
    indices_v = np.zeros((num_elements, 3), np.int_)
    i_local = np.zeros((num_elements, 3))
    v_local = np.zeros((num_elements, 3))
    for k in prange(num_elements):
        element = elements[k]
        vec_elem = value[k] / 12 * np.pi * determinant_b[element] * (bmat_11[element] * np.array([1, 2, 1]) +
                                                                     bmat_12[element] * np.array([1, 1, 2]) +
                                                                     r1[element] * np.array([4, 4, 4]))
        idx_node = elem2node[element, :]

        indices_v[k] = np.arange(3 * element, 3 * element + 3, 1)
        i_local[k] = idx_node
        v_local[k] = vec_elem

    # Write results to v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_load_elems_function(elements: np.ndarray, i: np.ndarray, v: np.ndarray,
                             weights: np.ndarray, local_coordinates: np.ndarray,
                             evaluation_points: np.ndarray, evaluations: np.ndarray,
                             determinant_b: np.ndarray,
                             elem2node) -> NoReturn:
    """Calculates elements for the load vector with load being a function.

    Calculates the contribution of the elements specified in `elements` to the load vector. The calculated values and
    line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in `elements`
        T       Number of elements in the mesh
        N       Number of evaluation points
        K       Size of i, j, v. Equal to 9T
        ======  =======

    Parameters
    ----------
    elements :  np.ndarray
        (E,) array. The indices of elements that should be considered.
    i : np.ndarray
        (K,) array. See the variable `i_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See the variable  `v_elem` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_element`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_element`.
    determinant_b : np.ndarray
        (T,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    elem2node : np.ndarray
        (T,3) array. See the filed in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_elements = len(elements)

    indices_v = np.zeros((num_elements, 3), np.int_)
    i_local = np.zeros((num_elements, 3))
    v_local = np.zeros((num_elements, 3))
    for k in prange(num_elements):
        element = elements[k]

        integral = np.zeros(3)
        for kk, weight in enumerate(weights):
            m_values = np.array([1 - local_coordinates[kk, 0] - local_coordinates[kk, 1], local_coordinates[kk, 0],
                                 local_coordinates[kk, 1]])
            integral = integral + weight * evaluations[k, kk] * evaluation_points[k, kk, 0] * m_values

        vec_elem = determinant_b[element] * 2 * np.pi * integral
        idx_node = elem2node[element, :]

        indices_v[k] = np.arange(3 * element, 3 * element + 3, 1)
        i_local[k] = idx_node
        v_local[k] = vec_elem

    # Write results to v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


# endregion


@njit_cache(parallel=True)
def calc_curl(val2edge: np.ndarray, node: np.ndarray, elem2node: np.ndarray, num_elem: int, b: np.ndarray,
              c: np.ndarray, s: np.ndarray):
    """Performant implementation of curl.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in the mesh
        N       Number of nodes in the mesh
        ======  =======

    Parameters
    ----------
    val2edge : np.ndarray
        (N,) array. One value per edge, the edge functions are defined on, i.e. one value per node in the mesh.
    node : np.ndarray
        (N,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    elem2node : np.ndarray
        (E,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    num_elem : int
        Total number of elements. See field in :py:class:`source.mesh.AxiMesh`.
    b : np.ndarray
        (E,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (E,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    s : s : np.ndarray
        (E,) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.

    Returns
    -------
    None
    """
    out = np.empty((num_elem, 2), dtype=val2edge.dtype)
    for k in prange(num_elem):
        nodes = elem2node[k]
        radius_mid = 1 / 3 * np.sum(node[nodes, 0])
        tmp = np.array([-1 / (2 * np.pi * s[k] * radius_mid) * np.sum(c[k] * val2edge[nodes]),
                        1 / (np.pi * s[k]) * np.sum(b[k] * val2edge[nodes])])
        out[k] = tmp
    return out


@njit_cache(parallel=True)
def calc_neumann_term_constant(edges: np.ndarray, i: np.ndarray, v: np.ndarray, value: Union[int, float],
                               node: np.ndarray, edge2node: np.ndarray):
    """Calculates elements for the neumann term vector with a constant value.

    Calculates the contribution of the elements specified in `edges` to the neumann term vector. The calculated values
    and line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 2*T
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the boundary condition is applied to.
    i : np.ndarray
        (K,) array. See variable `i_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See variable `v_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    value : Union[int, float]
        The value of the boundary condition.
    node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_edges = len(edges)
    indices_v = np.zeros((num_edges, 2), np.int_)
    i_local = np.zeros((num_edges, 2))
    v_local = np.zeros((num_edges, 2))
    for k in prange(num_edges):
        edge = edges[k]
        idx_node = edge2node[edge, :]
        coord_node1 = node[idx_node[0]]
        coord_node2 = node[idx_node[1]]
        length_edge = np.linalg.norm(coord_node2 - coord_node1)

        indices_v[k] = np.arange(2 * edge, 2 * edge + 2, 1)
        i_local[k] = idx_node
        v_local[k] = value * 2 * np.pi * length_edge * \
                     np.array([coord_node1[0] / 3 + coord_node2[0] / 6, coord_node1[0] / 6 + coord_node2[0] / 3])

    # Write results to i and v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_neumann_term_array(edges: np.ndarray, i: np.ndarray, v: np.ndarray, value: np.ndarray,
                            node: np.ndarray, edge2node: np.ndarray):
    """Calculates elements for the neumann term vector with a constant value.

    Calculates the contribution of the elements specified in `edges` to the neumann term vector. The calculated values
    and line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 2*T
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the boundary condition is applied to.
    i : np.ndarray
        (K,) array. See variable `i_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See variable `v_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    value : np.ndarray
        The values of the boundary condition.
    node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_edges = len(edges)

    indices_v = np.zeros((num_edges, 2), np.int_)
    i_local = np.zeros((num_edges, 2))
    v_local = np.zeros((num_edges, 2))
    for k in prange(num_edges):
        edge = edges[k]
        idx_node = edge2node[edge, :]
        coord_node1 = node[idx_node[0]]
        coord_node2 = node[idx_node[1]]
        length_edge = np.linalg.norm(coord_node2 - coord_node1)

        indices_v[k] = np.arange(2 * edge, 2 * edge + 2, 1)
        i_local[k] = idx_node
        v_local[k] = value[k] * 2 * np.pi * length_edge * \
                     np.array([coord_node1[0] / 3 + coord_node2[0] / 6, coord_node1[0] / 6 + coord_node2[0] / 3])

    # Write results to i and v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_neumann_term_function(edges: np.ndarray, i: np.ndarray, v: np.ndarray, weights: np.ndarray,
                               local_coordinates: np.ndarray, evaluation_points: np.ndarray,
                               evaluations: np.ndarray, node: np.ndarray,
                               edge2node: np.ndarray):
    """Calculates elements for the neumann term vector with a space dependant boundary condition.

    Calculates the contribution of the elements specified in `edges` to the neumann term vector. The calculated values
    and line information are written to v and i, respectively.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 2*T
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the face current is on.
    i : np.ndarray
        (K,) array. See variable `i_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    v : np.ndarray
        (K,) array. See variable `v_face` in :py:func:`TriAxisymmetricNodalShapeFunction.load_vector`.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all elements. See :py:func:`calc_evaluation_points_edge`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_edge`.
    node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    None
    """
    num_edges = len(edges)

    indices_v = np.zeros((num_edges, 2), np.int_)
    i_local = np.zeros((num_edges, 2))
    v_local = np.zeros((num_edges, 2))
    for k in range(num_edges):
        edge = edges[k]
        length_edge = np.linalg.norm(node[edge2node[edge, 1]] - node[edge2node[edge, 0]])
        integral = np.zeros(2)
        for kk, weight in enumerate(weights):
            m_values = np.array([1 - local_coordinates[kk], local_coordinates[kk]])
            integral = integral + weight * evaluations[k, kk] * 2 * np.pi * \
                       evaluation_points[k, kk, 0] * m_values * length_edge

        idx_node = edge2node[edge, :]

        indices_v[k] = np.arange(2 * edge, 2 * edge + 2, 1)
        i_local[k] = idx_node
        v_local[k] = integral

    # Write results to i and v
    idx = indices_v.flatten()
    i[idx] = i_local.flatten()
    v[idx] = v_local.flatten()


@njit_cache(parallel=True)
def calc_robin_term_constant(edges: np.ndarray, value: Union[int, float],
                             node: np.ndarray, edge2node: np.ndarray) -> Tuple[np.ndarray,
                                                                               np.ndarray,
                                                                               np.ndarray]:
    """Calculates elements for the robin term matrix with a constant boundary condition.

    Calculates the contribution of the elements specified in `edges` to the robin term matrix. the calculated values are
    returned.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 4*E
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the face current is on.
    value : Union[int, float]
        coefficient_dirichlet / coefficient_neumann
    node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    lines : np.ndarray
        (K,) array. Lines of the matrix elements
    columns : np.ndarray
        (K,) array. Columns of the matrix elements
    values : np.ndarray
        (K,) array. Values of the matrix elements
    """
    num_edges = len(edges)

    i_local = np.zeros((num_edges, 4))
    j_local = np.zeros((num_edges, 4))
    v_local = np.zeros((num_edges, 4))
    for k in prange(num_edges):
        edge = edges[k]
        idx_node = edge2node[edge, :]
        coord_node1 = node[idx_node[0]]
        coord_node2 = node[idx_node[1]]
        length_edge = np.linalg.norm(coord_node2 - coord_node1)

        mat_elem = 1 / 6 * np.pi * length_edge * value * (coord_node1[0] * np.array([[3, 1], [1, 1]]) +
                                                          coord_node2[0] * np.array([[1, 1], [1, 3]]))

        idx_node = edge2node[edge, :]
        tmp = np.stack((idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node))

        i_local[k] = np.reshape(tmp, (4,))
        j_local[k] = np.reshape(tmpt, (4,))
        v_local[k] = np.reshape(mat_elem, (4,))

    return i_local.flatten(), j_local.flatten(), v_local.flatten()


@njit_cache(parallel=True)
def calc_robin_term_array(edges: np.ndarray, value: np.ndarray,
                          node: np.ndarray, edge2node: np.ndarray) -> Tuple[np.ndarray,
                                                                            np.ndarray,
                                                                            np.ndarray]:
    """Calculates elements for the robin term matrix with a constant boundary condition.

    Calculates the contribution of the elements specified in `edges` to the robin term matrix. the calculated values are
    returned.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 4*E
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the face current is on.
    value : np.ndarray
        coefficient_dirichlet / coefficient_neumann
    node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    lines : np.ndarray
        (K,) array. Lines of the matrix elements
    columns : np.ndarray
        (K,) array. Columns of the matrix elements
    values : np.ndarray
        (K,) array. Values of the matrix elements
    """
    num_edges = len(edges)

    i_local = np.zeros((num_edges, 4))
    j_local = np.zeros((num_edges, 4))
    v_local = np.zeros((num_edges, 4))
    for k in prange(num_edges):
        edge = edges[k]
        idx_node = edge2node[edge, :]
        coord_node1 = node[idx_node[0]]
        coord_node2 = node[idx_node[1]]
        length_edge = np.linalg.norm(coord_node2 - coord_node1)

        mat_elem = 1 / 6 * np.pi * length_edge * value[k] * (coord_node1[0] * np.array([[3, 1], [1, 1]]) +
                                                             coord_node2[0] * np.array([[1, 1], [1, 3]]))

        idx_node = edge2node[edge, :]
        tmp = np.stack((idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node))

        i_local[k] = np.reshape(tmp, (4,))
        j_local[k] = np.reshape(tmpt, (4,))
        v_local[k] = np.reshape(mat_elem, (4,))

    return i_local.flatten(), j_local.flatten(), v_local.flatten()


@njit_cache(parallel=True)
def calc_robin_term_function(edges: np.ndarray, weights: np.ndarray, local_coordinates: np.ndarray,
                             evaluation_points: np.ndarray, evaluations: np.ndarray,
                             node_transformed: np.ndarray, edge2node: np.ndarray) -> Tuple[np.ndarray,
                                                                                           np.ndarray,
                                                                                           np.ndarray]:
    """Calculates elements for the robin term matrix with a space dependent boundary condition.

    Calculates the contribution of the elements specified in `edges` to the robin term matrix. the calculated values are
    returned.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of edges in `edges`
        T       Number of edges in the mesh
        N       Number of evaluation points
        K       Size of i, v. Equal to 4*E
        ======  =======

    Parameters
    ----------
    edges : np.ndarray
        (E,) array. Indices of the edges the face current is on.
    weights : np.ndarray
        (N,) array. The weights of the numerical integration.
    local_coordinates : np.ndarray
        (N,2) array. Two coordinates for each evaluation point.
    evaluation_points : np.ndarray
        (E,N,2) array. The evaluation points for all edges. See :py:func:`calc_evaluation_points_element`.
    evaluations : np.ndarray
        (E,N) array. For every element the N evalutations of the material function at the evaluation points for the
        numerical integration. See :py:func:`evaluate_for_num_int_element`.
    node_transformed : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    edge2node : np.ndarray
        (T,2) array. See field in :py:class:`source.mesh.AxiMesh`.

    Returns
    -------
    lines : np.ndarray
        (K,) array. Lines of the matrix elements
    columns : np.ndarray
        (K,) array. Columns of the matrix elements
    values : np.ndarray
        (K,) array. Values of the matrix elements
    """
    num_edges = len(edges)

    i_local = np.zeros((num_edges, 4))
    j_local = np.zeros((num_edges, 4))
    v_local = np.zeros((num_edges, 4))
    for k in prange(num_edges):
        edge = edges[k]

        coord_node1 = node_transformed[edge2node[edge, 0]]
        coord_node2 = node_transformed[edge2node[edge, 1]]
        length_edge = np.linalg.norm(coord_node2 - coord_node1)

        integral = np.zeros((2, 2))
        for kk, weight in enumerate(weights):
            m_values = np.array([1 - local_coordinates[kk], local_coordinates[kk]])
            integral = integral + weight * evaluations[k, kk] * \
                       evaluation_points[k, kk, 0] * np.outer(m_values, m_values) * length_edge

        mat_elem = 2 * np.pi * integral
        idx_node = edge2node[edge, :]

        tmp = np.stack((idx_node, idx_node))
        tmpt = np.column_stack((idx_node, idx_node))

        i_local[k] = np.reshape(tmp, (4,))
        j_local[k] = np.reshape(tmpt, (4,))
        v_local[k] = np.reshape(mat_elem, (4,))

    return i_local.flatten(), j_local.flatten(), v_local.flatten()


@njit_cache(parallel=True)
def calc_gradient(val2node: np.ndarray, elem2node: np.ndarray, num_elem: int, b: np.ndarray,
                  c: np.ndarray, element_areas_double: np.ndarray) -> np.ndarray:
    r"""Performant implementation of gradient.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        E       Number of elements in the mesh
        N       Number of nodes in the mesh
        ======  =======

    Parameters
    ----------
    val2node : np.ndarray
        (N,) array. One value per edge, the edge functions are defined on, i.e. one value per node in the mesh.
    elem2node : np.ndarray
        (E,2) array. See field in :py:class:`source.mesh.AxiMesh`.
    num_elem : int
        Total number of elements. See field in :py:class:`source.mesh.AxiMesh`.
    b : np.ndarray
        (E,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    c : np.ndarray
        (E,2) array. See the field in :py:class:`TriAxisymmetricNodalShapeFunction`.
    element_areas_double : s : np.ndarray
        (E,) array. Vector with 2 times the oriented area of all elements. See the field __parallelogram_areas_oriented
        in :py:class:`TriAxisymmetricNodalShapeFunction`.

    Returns
    -------
    None
    """
    out = np.empty((num_elem, 2), dtype=val2node.dtype)
    for k in prange(num_elem):
        area_double = element_areas_double[k]
        nodes = elem2node[k]
        tmp = np.array([np.sum(b[k] * val2node[nodes]),
                        np.sum(c[k] * val2node[nodes])])
        out[k] = tmp / area_double
    return out


class TriAxisymmetricNodalShapeFunction(NodalShapeFunction):
    """Class representing shape functions in axisymmetric geometries with a triangular mesh in the rz-plane."""

    def __init__(self, mesh: 'AxiMesh'):
        """Constructor.

        Parameters
        ----------
        mesh : AxiMesh
            The mesh object.
        """
        super().__init__(mesh, dim=2, allocation_size=9 * mesh.num_elem)
        self.__determinant_b = np.zeros(mesh.num_elem)
        self.__calc_coefficients()

    def _calc_matrix_constant_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: float, weights: np.ndarray, local_coordinates: np.ndarray,
                                     *args, evaluator: str = None, **kwargs):
        if matrix_type == "divgrad":
            calc_divgrad_constant_scalar(indices, i, j, v, value, self.b, self.c, self.bmat_11, self.bmat_12, self.r1,
                                         self.__determinant_b, self.mesh.elem2node)
        elif matrix_type == "mass":
            calc_mass_constant_scalar(indices, i, j, v, value, self.bmat_11, self.bmat_12, self.r1,
                                      self.__determinant_b, self.mesh.elem2node)
        else:
            raise NotImplementedError(f"calc_matrix_constant_scalar is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_constant_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        if matrix_type == "divgrad":
            calc_divgrad_constant_tensor(indices, i, j, v, value, self.b, self.c, self.bmat_11, self.bmat_12, self.r1,
                                         self.__determinant_b, self.mesh.elem2node)
        elif matrix_type == "mass":
            raise NotImplementedError('Mass matrix does not support tensors.')
        else:
            raise NotImplementedError(f"calc_matrix_constant_tensor is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_scalar_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        if evaluator == "eval_hom_nonlin_iso":
            evaluations = eval_hom_nonlin_iso(value, indices)
        elif evaluator is None:
            evaluations = value.astype(float)
        else:
            raise NotImplementedError(f"Evaluator {evaluator} is not implemented.")

        if matrix_type == "divgrad":
            calc_divgrad_scalar_per_elem(indices, i, j, v, evaluations, self.b, self.c,
                                         self.bmat_11, self.bmat_12, self.r1, self.__determinant_b, self.mesh.elem2node)
        elif matrix_type == "mass":
            calc_mass_scalar_per_elem(indices, i, j, v, evaluations, self.bmat_11, self.bmat_12,
                                      self.r1, self.__determinant_b, self.mesh.elem2node)
        else:
            raise NotImplementedError(f"calc_matrix_scalar_per_elem is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_tensor_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray, weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        if evaluator == "eval_hom_nonlin_aniso":
            evaluations = eval_hom_nonlin_aniso(value, indices)
        elif evaluator is None:
            evaluations = value.astype(float)
        else:
            raise NotImplementedError(f"Evaluator {evaluator} is not implemented.")

        if matrix_type == "divgrad":
            calc_divgrad_tensor_per_elem(indices, i, j, v, evaluations, self.b, self.c,
                                         self.bmat_11, self.bmat_12, self.r1, self.__determinant_b, self.mesh.elem2node)
        elif matrix_type == "mass":
            raise NotImplementedError('Mass matrix does not support tensors.')
        else:
            raise NotImplementedError(f"calc_matrix_tensor_per_elem is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_function_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., float], weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        evaluation_points = calc_evaluation_points_element(local_coordinates, self.transform_coefficients,
                                                           indices)
        if evaluator == "eval_inhom_lin_iso":
            evaluations = eval_inhom_lin_iso(value, evaluation_points)
        elif evaluator == "eval_inhom_nonlin_iso":
            evaluations = eval_inhom_nonlin_iso(value, evaluation_points, indices)
        else:
            raise NotImplementedError(f"Evaluator {evaluator} is not implemented.")

        if matrix_type == "divgrad":
            calc_divgrad_function_scalar(indices, i, j, v, weights, evaluation_points, evaluations,
                                         self.b, self.c, self.__determinant_b,
                                         self.mesh.elem2node)
        elif matrix_type == "mass":
            calc_mass_function_scalar(indices, i, j, v, weights, local_coordinates, evaluation_points, evaluations,
                                      self.__determinant_b, self.mesh.elem2node)
        else:
            raise NotImplementedError(f"calc_matrix_function_scalar is not implemented for matrix type: {matrix_type}.")

    def _calc_matrix_function_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., np.ndarray], weights: np.ndarray,
                                     local_coordinates: np.ndarray, *args, evaluator: str = None, **kwargs) -> None:
        evaluation_points = calc_evaluation_points_element(local_coordinates, self.transform_coefficients,
                                                           indices)
        if evaluator == "eval_inhom_lin_aniso":
            evaluations = eval_inhom_lin_aniso(value, evaluation_points)
        elif evaluator == "eval_inhom_nonlin_aniso":
            evaluations = eval_inhom_nonlin_aniso(value, evaluation_points, indices)
        else:
            raise NotImplementedError(f"Evaluator {evaluator} is not implemented.")

        if matrix_type == "divgrad":
            calc_divgrad_function_tensor(indices, i, j, v, weights, evaluation_points, evaluations,
                                         self.b, self.c, self.__determinant_b, self.mesh.elem2node)
        elif matrix_type == "mass":
            raise NotImplementedError('Mass matrix does not support tensors.')
        else:
            raise NotImplementedError(f"calc_matrix_function_scalar is not implemented for matrix type: {matrix_type}.")

    # @staticmethod
    # def base1(r: float, z: float) -> float:
    #     """
    #     First basis function on reference element for the mass matrix.
    #
    #     Parameters
    #     ----------
    #     r : float
    #         r coordinate
    #     z : float
    #         z coordinate
    #
    #     Returns
    #     -------
    #     float
    #
    #     Notes
    #     -----
    #     Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf
    #     """
    #     # return 1 / 3 * (4 - r ** 2 - 3 * z)
    #     return 1 - r - z
    #
    # # noinspection PyUnusedLocal
    # @staticmethod
    # def base2(r: float, z: float) -> float:
    #     """
    #     Second basis function on reference element for the mass matrix.
    #
    #     Parameters
    #     ----------
    #     r : float
    #         r coordinate
    #     z : float
    #         z coordinate
    #
    #     Returns
    #     -------
    #     float
    #
    #     Notes
    #     -----
    #     Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf
    #     """
    #     # return 1 / 3 * (-1 + r ** 2)
    #     return r
    #
    # # noinspection PyUnusedLocal
    # @staticmethod
    # def base3(r: float, z: float) -> float:
    #     """
    #     Third basis function on reference element for the mass matrix.
    #
    #     Parameters
    #     ----------
    #     r : float
    #         r coordinate
    #     z : float
    #         z coordinate
    #
    #     Returns
    #     -------
    #     float
    #
    #     Notes
    #     -----
    #     Explained in info/Ansatzfunktionen/ansatzfunktionen.pdf
    #     """
    #     return z

    def __calc_coefficients(self) -> NoReturn:
        """
        Calculates some coefficients that are needed in various methods of this class.

        Returns
        -------
        NoReturn
        """
        r1 = self.mesh.node[self.mesh.elem2node, 0]
        z1 = self.mesh.node[self.mesh.elem2node, 1]
        r2 = np.roll(r1, -1, axis=1)
        z2 = np.roll(z1, -1, axis=1)
        r3 = np.roll(r1, -2, axis=1)
        z3 = np.roll(z1, -2, axis=1)
        self.a = r2 * z3 - r3 * z2
        self.b = z2 - z3
        self.c = r3 - r2

        r1 = self.mesh.node[self.mesh.elem2node[:, 0], 0]
        r2 = self.mesh.node[self.mesh.elem2node[:, 1], 0]
        r3 = self.mesh.node[self.mesh.elem2node[:, 2], 0]
        z1 = self.mesh.node[self.mesh.elem2node[:, 0], 1]
        z2 = self.mesh.node[self.mesh.elem2node[:, 1], 1]
        z3 = self.mesh.node[self.mesh.elem2node[:, 2], 1]
        self.p1_to_p2 = np.c_[r2 - r1, z2 - z1]
        self.p1_to_p3 = np.c_[r3 - r1, z3 - z1]
        # (elements,2,3) array. For every element there are the 3 coefficients for s and z
        self.transform_coefficients = np.stack((self.p1_to_p2, self.p1_to_p3, np.c_[r1, z1]), axis=2)

        self.bmat_11 = self.p1_to_p2[:, 0]
        self.bmat_12 = self.p1_to_p3[:, 0]
        self.r1 = r1

        # Calculate the cross product of the vectors p1_to_p2 and p1_to_p3. The absolute value corresponds to the area
        # of the parallelogram spanned by the two vector (= twice the element area) and the sign depends on the
        # orientation of the triangle (ijk or ikj)
        self.__parallelogram_areas_oriented = self.p1_to_p2[:, 0] * self.p1_to_p3[:, 1] - \
            self.p1_to_p2[:, 1] * self.p1_to_p3[:, 0]
        # Calculate the absolute value of the determinant of the transformation matrix for every element
        self.__determinant_b = np.abs(self.__parallelogram_areas_oriented)

    def divgrad_operator(self, *material: Union[Callable[[float, float], float],
                                                ndarray, float, Union['Regions', 'Materials', 'MatProperty']],
                         integration_order: int = 3) -> coo_matrix:
        return self._matrix_routine("divgrad", *material, integration_order=integration_order)

    def mass_matrix(self, *material: Union[Callable[[float, float], float],
                                           ndarray, float, Union['Regions', 'Materials', 'MatProperty']],
                    integration_order: int = 3) -> coo_matrix:
        return self._matrix_routine("mass", *material, integration_order=integration_order)

    def load_vector(self, *load: Union[Callable[[float, float], float], ndarray, float,
                                       Union['Regions', 'Excitations']], integration_order: int = 3) -> coo_matrix:
        case, load = self._process_load(*load)

        if not isinstance(integration_order, int):
            raise ValueError(f"integration order is of type {type(integration_order)} but has to be an int.")

        weights, local_coordinates = gauss_triangle(integration_order)
        weights_edge, local_coordinates_edge = gauss(integration_order)

        # For i_node and v_node a list is taken because the number of line currents is most likely low compared to the
        # number of all nodes
        i_node = []
        v_node = []
        i_face = np.zeros(2 * self.mesh.num_edge, dtype=np.int32)
        v_face = np.zeros(2 * self.mesh.num_edge)
        i_elem = np.zeros(3 * self.mesh.num_elem, dtype=np.int32)
        v_elem = np.zeros(3 * self.mesh.num_elem)

        flag_node = True
        flag_face = False

        if case == "tuple":  # material is a tuple
            regions, excitations = load
            for key in regions.get_keys():
                regi = regions.get_regi(key)

                exci_id = regi.exci
                if exci_id is None:
                    continue
                exci = excitations.get_exci(exci_id)

                if regi.dim == 0:
                    flag_node = True
                    indices = np.where(self.mesh.node2regi == regi.ID)[0]
                    for k in indices:
                        i_node.append(k)
                        if exci.is_constant:
                            v_node.append(exci.value * 2 * np.pi * self.mesh.node[k, 0])
                        elif not exci.is_homogeneous:
                            v_node.append(exci.value(self.mesh.node[k, 0], self.mesh.node[k, 1]) *
                                          2 * np.pi * self.mesh.node[k, 0])
                        elif exci.is_time_dependent or not exci.is_linear:
                            v_node.append(exci.value() * 2 * np.pi * self.mesh.node[k, 0])
                        else:  # When the point source depends on time or any field
                            raise NotImplementedError
                elif regi.dim == 1:
                    flag_face = True

                    indices = np.where(self.mesh.edge2regi == regi.ID)[0]
                    try:
                        indices.astype(np.int32, casting="safe", copy=False)
                    except TypeError:
                        pass

                    if exci.is_constant:
                        calc_neumann_term_constant(indices, i_face, v_face, exci.value, self.mesh.node,
                                                   self.mesh.edge2node)
                    elif not exci.is_homogeneous:
                        evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node,
                                                                        indices, self.mesh.edge2node)
                        evaluations = evaluate_for_num_int_edge(exci.value, evaluation_points)
                        calc_neumann_term_function(indices, i_face, v_face, weights_edge, local_coordinates_edge,
                                                   evaluation_points, evaluations, self.mesh.node,
                                                   self.mesh.edge2node)
                    elif not exci.is_linear or exci.is_time_dependent:
                        value = exci.value()
                        if len(value) == 1:
                            calc_neumann_term_constant(indices, i_face, v_face, value, self.mesh.node,
                                                       self.mesh.edge2node)
                        elif len(value) == len(indices):
                            calc_neumann_term_array(indices, i_face, v_face, value, self.mesh.node,
                                                    self.mesh.edge2node)
                        else:
                            raise NotImplementedError

                    else:  # When the line source depends on time or any field
                        raise NotImplementedError
                elif regi.dim == 2:
                    indices = np.where(self.mesh.elem2regi == regi.ID)[0]
                    try:
                        indices.astype(np.int32, casting="safe", copy=False)
                    except TypeError:
                        pass

                    if exci.is_constant:
                        calc_load_elems_constant(indices, i_elem, v_elem, exci.value,
                                                 self.bmat_11, self.bmat_12, self.r1, self.__determinant_b,
                                                 self.mesh.elem2node)
                    elif not exci.is_constant:
                        if not exci.is_homogeneous:
                            evaluation_points = calc_evaluation_points_element(local_coordinates,
                                                                               self.transform_coefficients, indices)
                            evaluations = eval_inhom_lin_iso(exci.value, evaluation_points)
                            calc_load_elems_function(indices, i_elem, v_elem, weights, local_coordinates,
                                                     evaluation_points,
                                                     evaluations, self.__determinant_b, self.mesh.elem2node)

                        elif exci.is_time_dependent or (not exci.is_linear):  # change: hacky
                            calc_load_elems_constant(indices, i_elem, v_elem, exci.value(), weights, local_coordinates,
                                                     self.transform_coefficients, self.__determinant_b,
                                                     self.mesh.elem2node)

                    else:
                        raise NotImplementedError
        elif case == "number":
            calc_load_elems_constant(np.arange(self.mesh.num_elem), i_elem, v_elem, load,
                                     self.bmat_11, self.bmat_12, self.r1, self.__determinant_b, self.mesh.elem2node)
        elif case == "array":  # material is an array
            calc_load_elems_vector(np.arange(self.mesh.num_elem), i_elem, v_elem, load,
                                   self.bmat_11, self.bmat_12, self.r1, self.__determinant_b, self.mesh.elem2node)
        elif case == "function":  # material is a function
            evaluation_points = calc_evaluation_points_element(local_coordinates, self.transform_coefficients,
                                                               np.arange(self.mesh.num_elem))
            evaluations = eval_inhom_lin_iso(load, evaluation_points)
            calc_load_elems_function(np.arange(self.mesh.num_elem), i_elem, v_elem, weights, local_coordinates,
                                     evaluation_points, evaluations,
                                     self.__determinant_b, self.mesh.elem2node)
        else:
            raise ValueError("Argument type not expected. Function can not handle the case of a " + case)
        i = i_elem
        v = v_elem
        if flag_node:
            i = np.concatenate((i, np.array(i_node, dtype=int)))
            v = np.concatenate((v, np.array(v_node)))
        if flag_face:
            i = np.concatenate((i, i_face))
            v = np.concatenate((v, v_face))

        columns = np.zeros(len(i), dtype=int)
        load_vector = coo_matrix((v, (i, columns)), shape=(self.mesh.num_node, 1))
        return load_vector

    # pylint: disable=line-too-long
    def neumann_term(self, *args: Union[Tuple['Regions', 'BdryCond'],
                                        Tuple[ndarray, Union[Tuple[Union[float, ndarray, Callable[..., float]], ...],
                                                             Callable[..., Tuple[float]]]]],
                     integration_order: int = 3) -> coo_matrix:
        r"""Compute the Neumann term on a part of the boundary (see the notes).

        Parameters
        ----------
        args : Union[Tuple[Regions, BdryCond], Tuple[ndarray, Union[Tuple[Union[float, ndarray, Callable[...,
        float]], ...], Callable[..., Tuple[float]]]]]
            The following possibilities are:

            - **Regions, BdyCond**: Will search for instances of BCNeumann in BdryCond and use the value of these.
            - **Tuple[ndarray, Union[Tuple[Union[float, ndarray, Callable[..., float]], ...], Callable[..., Tuple[float]]]]**:
              The first argument is an array of shape (N,) and contains the indices of the boundary elements. The second
              argument is the value of g at these elements.
        integration_order : int, optional
            The integration order for the case that a numerical integration is necessary. Default is 1

        Returns
        -------
        q : csr_matrix
            The vector :math:`\mathbf{q}` for the Neumann term.

        Notes
        -----
        The Neumann term for edge basis functions is a vector :math:`\mathbf{q}` where the elements are defined as

        .. math::

            q_i = \int_{\partial\Omega} \vec{w}_i\cdot\vec{g}\,\mathrm{d} S\,.

        In the finite element formulation it is :math:`\vec{g} = \vec{n}\times\vec{H}`.
        """
        flag_regions, flag_value = self._process_neumann(*args, allow_indices_tuple=True)

        weights, local_coordinates = gauss(integration_order)  # quadrature along unit line
        i = np.zeros(2 * self.mesh.num_edge, dtype=np.int32)
        v = np.zeros(2 * self.mesh.num_edge)

        if flag_regions:  # We have regions and boundary_conditions
            regions: Regions = args[0]
            boundary_conditions: BdryCond = args[1]
            for key in boundary_conditions.get_ids():
                bc = boundary_conditions.get_bc(key)
                if not isinstance(bc, BCNeumann):
                    continue

                indices_regions = regions.find_regions_of_boundary_condition(bc.ID)
                indices_edges = np.empty(0)  # Indices of all edges that are on the boundary bc
                for region_id in indices_regions:
                    indices_edges = np.r_[indices_edges, np.where(self.mesh.edge2regi == region_id)[0]]
                indices_edges = indices_edges.astype(int)

                if not bc.is_homogeneous:
                    evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node,
                                                                    indices_edges, self.mesh.edge2node)
                    evaluations = evaluate_for_num_int_edge(bc.value, evaluation_points)
                    calc_neumann_term_function(indices_edges, i, v, weights, local_coordinates, evaluation_points,
                                               evaluations, self.mesh.node, self.mesh.edge2node)

                elif bc.is_constant:
                    calc_neumann_term_constant(indices_edges, i, v, bc.value, self.mesh.node, self.mesh.edge2node)
                else:
                    calc_neumann_term_constant(indices_edges, i, v, bc.value(), self.mesh.node, self.mesh.edge2node)

        else:  # We have indices and value
            indices, value = args[0], args[1][0] if isinstance(args[1], Tuple) else args[1]

            if flag_value == 'callable':
                evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node, indices,
                                                                self.mesh.edge2node)
                evaluations = evaluate_for_num_int_edge(value, evaluation_points)
                calc_neumann_term_function(indices, i, v, weights, local_coordinates, evaluation_points, evaluations,
                                           self.mesh.node, self.mesh.edge2node)
            elif flag_value == 'array':
                calc_neumann_term_array(indices, i, v, value, self.mesh.node, self.mesh.edge2node)
            elif flag_value == 'value':
                calc_neumann_term_constant(indices, i, v, value, self.mesh.node, self.mesh.edge2node)
        columns = np.zeros(len(i), dtype=int)
        neumann_vector = coo_matrix((v, (i, columns)), shape=(self.mesh.num_node, 1))
        return neumann_vector

    # pylint: disable=cell-var-from-loop
    def robin_terms(self, regions: 'Regions', boundary_condition: 'BdryCond', integration_order: int = 1) -> \
            Tuple[csr_matrix, csr_matrix]:
        bc_keys = boundary_condition.get_ids()
        weights, local_coordinates = gauss(integration_order)
        values_list = []
        lines_list = []
        columns_list = []
        robin_vector = coo_matrix((self.mesh.num_node, 1))

        for bc_key in bc_keys:
            bc = boundary_condition.get_bc(bc_key)
            if not isinstance(bc, BCRobin):
                continue

            # Calculate the matrix
            indices_regions = regions.find_regions_of_boundary_condition(bc.ID)
            indices_edges = np.empty(0)  # Indices of all edges that are on the boundary bc
            for region_id in indices_regions:
                indices_edges = np.r_[indices_edges, np.where(self.mesh.edge2regi == region_id)[0]]
            indices_edges = indices_edges.astype('int')

            coef_dir = bc.coefficient_dirichlet
            coef_neum = bc.coefficient_neumann
            if isinstance(coef_neum, ndarray) or isinstance(coef_dir, ndarray):
                raise NotImplementedError
            value_bc = bc.value
            if callable(coef_dir):
                if callable(coef_neum):  # Both coefficients are functions
                    evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node,
                                                                    indices_edges, self.mesh.edge2node)
                    evaluations = evaluate_for_num_int_edge(lambda r, z: coef_dir(r, z) / coef_neum(r, z),
                                                            evaluation_points)
                    robin_lines_per_bc, robin_columns_per_bc, robin_values_per_bc = calc_robin_term_function(
                        indices_edges, weights, local_coordinates, evaluation_points, evaluations,
                        self.mesh.node, self.mesh.edge2node)
                else:  # Only dirichlet coefficient is a function
                    evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node,
                                                                    indices_edges, self.mesh.edge2node)
                    evaluations = evaluate_for_num_int_edge(lambda r, z: coef_dir(r, z) / coef_neum,
                                                            evaluation_points)
                    robin_lines_per_bc, robin_columns_per_bc, robin_values_per_bc = calc_robin_term_function(
                        indices_edges, weights, local_coordinates, evaluation_points, evaluations,
                        self.mesh.node, self.mesh.edge2node)
            else:
                if callable(coef_neum):  # Only coef_neum is a function
                    evaluation_points = calc_evaluation_points_edge(local_coordinates, self.mesh.node,
                                                                    indices_edges, self.mesh.edge2node)
                    evaluations = evaluate_for_num_int_edge(lambda r, z: coef_dir / coef_neum(r, z), evaluation_points)
                    robin_lines_per_bc, robin_columns_per_bc, robin_values_per_bc = calc_robin_term_function(
                        indices_edges, weights, local_coordinates, evaluation_points, evaluations,
                        self.mesh.node, self.mesh.edge2node)
                else:  # Both are no functions
                    value = coef_dir / coef_neum
                    robin_lines_per_bc, robin_columns_per_bc, robin_values_per_bc = calc_robin_term_constant(
                        indices_edges, value, self.mesh.node, self.mesh.edge2node)
            values_list.append(robin_values_per_bc)
            lines_list.append(robin_lines_per_bc)
            columns_list.append(robin_columns_per_bc)

            # Calculate the vector
            if callable(value_bc):
                if callable(coef_neum):
                    # Both are functions
                    robin_vector = robin_vector + self.neumann_term(indices_edges,
                                                                    (lambda r, z: value_bc(r, z) / coef_neum(r, z),),
                                                                    integration_order=integration_order)
                else:
                    # only value_bc is a function
                    if isinstance(coef_neum, np.ndarray):
                        robin_vector = robin_vector + self.neumann_term(indices_edges,
                                                                        (lambda r, z: value_bc(r, z) / coef_neum,),
                                                                        integration_order=integration_order)
                    elif isinstance(coef_neum, (float, int)):
                        robin_vector = robin_vector + self.neumann_term(indices_edges,
                                                                        (lambda r, z: value_bc(r, z) / coef_neum,),
                                                                        integration_order=integration_order)
                    else:
                        raise Exception("Type of coef_neum is not supported.")
            else:
                if callable(coef_neum):
                    # Only Neumann coefficient is a function
                    if isinstance(value_bc, np.ndarray):
                        robin_vector = robin_vector + self.neumann_term(indices_edges,
                                                                        (lambda r, z: value_bc / coef_neum(r, z),),
                                                                        integration_order=integration_order)
                    elif isinstance(value_bc, (float, int)):
                        robin_vector = robin_vector + self.neumann_term(indices_edges,
                                                                        (lambda r, z: value_bc / coef_neum(r, z),),
                                                                        integration_order=integration_order)
                    else:
                        raise Exception("Type of value not supported")

                else:
                    value_tmp = value_bc / coef_neum
                    if isinstance(value_tmp, np.ndarray):
                        robin_vector = robin_vector + self.neumann_term(indices_edges, (value_tmp,),
                                                                        integration_order=integration_order)
                    elif isinstance(value_tmp, (float, int)):
                        robin_vector = robin_vector + self.neumann_term(indices_edges, (value_tmp,),
                                                                        integration_order=integration_order)
                    else:
                        raise Exception("Type not supported.")

        robin_matrix = coo_matrix(
            (np.concatenate(values_list), (np.concatenate(lines_list), np.concatenate(columns_list))),
            shape=(self.mesh.num_node, self.mesh.num_node))
        return robin_matrix, robin_vector.tocoo()

    def gradient(self, val2node: np.ndarray) -> np.ndarray:
        out = calc_gradient(val2node, self.mesh.elem2node, self.mesh.num_elem, self.b, self.c,
                            self.__parallelogram_areas_oriented)
        return out
