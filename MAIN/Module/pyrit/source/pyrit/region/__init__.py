# coding=utf-8
"""
============================
Region (:mod:`pyrit.region`)
============================

.. currentmodule:: pyrit.region

This module organizes classes for the organization of materials, boundary conditions and excitations on the mesh are
provided. Regions can be seen as discrete counterpart to physical groups and are defined on the mesh. The package
includes a class for defining a region and a container class that manages different regions.

.. autosummary::
    :toctree: autosummary/

    Regi
    Regions

"""

from .Regi import Regi
from .Regions import Regions

__all__ = ['Regi', 'Regions']
