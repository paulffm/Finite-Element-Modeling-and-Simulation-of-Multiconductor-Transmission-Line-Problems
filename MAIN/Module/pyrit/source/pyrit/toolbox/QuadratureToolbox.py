# coding=utf-8
"""
Quadrature Toolbox

.. sectionauthor:: bundschuh, menzenbach
"""
# pylint: disable=unspecified-encoding

import pathlib
import json
from typing import Tuple
import numpy as np


file_directory = pathlib.Path(__file__).parent


def gauss(order: int, left: float = 0, right: float = 1) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Returns weights and points for Gauss-Legendre quadrature on the interval [left,right].

    Parameters
    ----------
    order : int
        The order of the quadrature rule.
    left : float, optional
        Left border of the interval. Default is 0
    right : float, optional
        Right border of the interval. Default is 1

    Returns
    -------
    weights : np.ndarray
        (N,) array. The weights of the quadrature rule.
    points : np.ndarray
        (N,) array. The points of the quadrature rule.

    Notes
    -----
    The weights :math:`w_i` and points :math:`x_i` are used as follows:
    :math:`\int_a^b f(x) dx \approx \sum_{i=1}{N} w_i f(x_i)`
    """
    if left >= right:
        raise ValueError("The left border has to be smaller than the right border.")

    try:
        p, w = np.polynomial.legendre.leggauss(order)
    except ValueError as ve:
        raise ValueError(f"Order is {order}, but has to be a positive integer.") from ve
    except TypeError as te:
        raise TypeError(f"Order is {order}, but has to be a positive integer.") from te

    weights = (right - left) / 2 * w
    points = (right - left) / 2 * p + (left + right) / 2

    return weights, points


def gauss_triangle(p: int) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Returns weights and local coordinates in reference triangle (spanned by (0,0), (0,1), (1,0)).

    Parameters
    ----------
    p : int
        Degree of the approximation.

    Returns
    -------
    weights : np.ndarray
        (N,) array with N weights
    local_cds : np.ndarray
        (N,2) array with N points in reference triangle and 2 coordinates (x and y) for each point

    Raises
    ------
    ValueError
        If the parameter p has the wrong type or is not available.

    Notes
    -----
    With T being the reference triangle, the weights :math:`w_i` and coordinates :math:`x_i` and :math:`y_i` are used as
    follows: :math:`\int_T f(x,y) dA \approx \sum_{i=1}^{N} w_i f(x_i,y_i)`.

    Points and weights from https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/rules_triangle.html
    """
    if not isinstance(p, int):
        raise ValueError("Order has to be an integer.")
    if not 1 <= p <= 20:
        raise ValueError(f"The order '{p}' is not available.")
    with open(file_directory.joinpath('schemes/dunavant_all.json'), 'r') as f:
        data = json.load(f)

    data = data[f'{p}']

    return np.atleast_1d(np.array(data['weights'])), np.atleast_2d(np.array(data['points']))


def gauss_tetrahedron(p: int) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Returns weights and local coordinates in reference tetrahedron (spanned by (0,0,0), (1,0,0), (0,1,0), (0,0,1)).

    Parameters
    ----------
    p : int
        Degree of the approximation.

    Returns
    -------
    weights : np.ndarray
        (N,) array with N weights
    points : np.ndarray
        (N,3) array with N points in reference tetrahedron and 3 coordinates (x, y and z) for each point

    Raises
    ------
    ValueError
        If the parameter p has the wrong type or is not available.

    Notes
    -----
    With T being the reference tetrahedron, the weights :math:`w_i` and coordinates :math:`x_i`, :math:`y_i` and
    :math:`z_i` are used as follows:

    .. :math:
        \int_T f(x, y, z) dV \approx \sum_{i=1}{N} w_i f(x_i, y_i, z_i)\,.

    Points and weights from https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6313
    """
    if not isinstance(p, int):
        raise ValueError("Order has to be an integer.")
    if not 1 <= p <= 20:
        raise ValueError(f"The order '{p}' is not available.")
    with open(file_directory.joinpath('schemes/jaskowiec_sukumar_all.json'), 'r') as f:
        data = json.load(f)

    data = data[f'{p}']

    return np.atleast_1d(np.array(data['weights'])), np.atleast_2d(np.array(data['points']))
