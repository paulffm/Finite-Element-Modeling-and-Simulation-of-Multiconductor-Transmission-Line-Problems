# coding=utf-8
"""
=============================================
Geo library (:mod:`pyrit.toolbox.geolibrary`)
=============================================

.. currentmodule:: pyrit.toolbox.geolibrary

This module contains a collection of :py:class:`Geo`-objects.

Primitive Geos
--------------

.. autosummary::
    :toctree: autosummary/

    FoilWinding
    Yoke


Complex Geos
------------

They consist of other Geos

.. autosummary::
    :toctree: autosummary/

    FoilWindingWithSurrounding
    PotInductor
"""
from .Yoke import Yoke
from .FoilWinding import FoilWinding
from .FoilWindingWithSurrounding import FoilWindingWithSurrounding
from .PotInductor import PotInductor

__all__ = ['FoilWinding', 'Yoke', 'FoilWindingWithSurrounding', 'PotInductor']
