# coding=utf-8
"""
================================
Material (:mod:`pyrit.material`)
================================

.. currentmodule:: pyrit.material

This module handles materials.

Materials and container class
-----------------------------

.. autosummary::
    :toctree: autosummary/

    Materials
    Conductivity
    Density
    VolumetricHeatCapacity
    Permeability
    Permittivity
    Reluctivity
    ThermalConductivity
    Resistivity
    DifferentialConductivity
    DifferentialPermittivity

Abstract classes
----------------

.. autosummary::
    :toctree: autosummary/

    Mat
    MatProperty

Exceptions
----------

.. autosummary::
    :toctree: autosummary/

    PropertyConversionException
"""

from .MatProperty import MatProperty
from .Conductivity import Conductivity
from .Density import Density
from .VolumetricHeatCapacity import VolumetricHeatCapacity
from .Permeability import Permeability
from .Permittivity import Permittivity
from .Reluctivity import Reluctivity
from .ThermalConductivity import ThermalConductivity
from .Resistivity import Resistivity
from .DifferentialConductivity import DifferentialConductivity
from .DifferentialPermittivity import DifferentialPermittivity
from .Mat import Mat, PropertyConversionException

from .Materials import Materials

__all__ = ['Conductivity', 'Density', 'VolumetricHeatCapacity', 'Mat', 'MatProperty', 'Materials',
           'Permeability', 'Permittivity', 'Reluctivity', 'ThermalConductivity', 'Resistivity',
           'DifferentialConductivity', 'DifferentialPermittivity', 'PropertyConversionException']
