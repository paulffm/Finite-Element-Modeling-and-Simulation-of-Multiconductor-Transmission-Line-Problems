# coding=utf-8
"""
===========================================
Shape function (:mod:`pyrit.shapefunction`)
===========================================

.. currentmodule:: pyrit.shapefunction

This package provides the core of the FE method, i.e. the different shape functions. There is one class per FE shape
function and dimension of the domain. For example, there is a class for nodal shape functions in axisymmetric
coordinates or a class for edge shape functions in two dimensional Cartesian coordinates. These classes implement the
routines to compute the FE matrices.

Actual shape function classes
-----------------------------

.. autosummary::
    :toctree: autosummary/

    TriCartesianNodalShapeFunction
    TriAxisymmetricNodalShapeFunction
    TriAxisymmetricEdgeShapeFunction
    TetCartesianNodalShapeFunction


Abstract class
--------------

.. autosummary::
    :toctree: autosummary/

    ShapeFunction
    NodalShapeFunction
    EdgeShapeFunction

"""

from .ShapeFunction import ShapeFunction
from .NodalShapeFunction import NodalShapeFunction
from .EdgeShapeFunction import EdgeShapeFunction
from .TriCartesianNodalShapeFunction import TriCartesianNodalShapeFunction
from .TriAxisymmetricNodalShapeFunction import TriAxisymmetricNodalShapeFunction
from .TriAxisymmetricEdgeShapeFunction import TriAxisymmetricEdgeShapeFunction
from .TetCartesianNodalShapeFunction import TetCartesianNodalShapeFunction
from .TriCartesianEdgeShapeFunction import TriCartesianEdgeShapeFunction

__all__ = ['ShapeFunction', 'NodalShapeFunction', 'EdgeShapeFunction', 'TriCartesianNodalShapeFunction',
           'TriAxisymmetricNodalShapeFunction', 'TriAxisymmetricEdgeShapeFunction', 'TetCartesianNodalShapeFunction',
           'TriCartesianEdgeShapeFunction']
