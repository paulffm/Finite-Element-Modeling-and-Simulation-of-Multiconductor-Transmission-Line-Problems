# coding=utf-8
"""
================================
Geometry (:mod:`pyrit.geometry`)
================================

.. currentmodule:: pyrit.geometry

The main part of the *Geometry* module is the class :py:class:`Geometry`. It is responsible for the connection between
**Pyrit** and *Gmsh*.

Main classes
------------

.. autosummary::
    :toctree: autosummary/

    Geometry
    PhysicalGroup

Classes for geometrical entities
--------------------------------

**0 Dimensional**

.. autosummary::
    :toctree: autosummary/

    Point


**1 Dimensional**

.. autosummary::
    :toctree: autosummary/

    Curve
    Line
    Circle
    CircleArc
    Ellipse
    EllipseArc

**2 Dimensional**

.. autosummary::
    :toctree: autosummary/

    Surface
    Rectangle

**3 Dimensional**

.. note::

    Not implemented yet.

Abstract classes
----------------

.. autosummary::
    :toctree: autosummary/

    PrimitiveGeo
    Geo
    BasicAnchor

"""

from .PrimitiveGeo import PrimitiveGeo
from .Point import Point
from .Curve import Curve
from .Surface import Surface

from .Line import Line
from .Circle import Circle
from .CircleArc import CircleArc
from .Ellipse import Ellipse
from .EllipseArc import EllipseArc

from .Rectangle import Rectangle

from .PhysicalGroup import PhysicalGroup
from .Geo import Geo, NotAddedToGmshError, BasicAnchor
from .Geometry import Geometry

__all__ = ['Circle', 'CircleArc', 'Curve', 'Ellipse', 'EllipseArc', 'Geo', 'Geometry', 'Line', 'PhysicalGroup', 'Point',
           'PrimitiveGeo', 'Rectangle', 'Surface', 'NotAddedToGmshError', 'BasicAnchor']
