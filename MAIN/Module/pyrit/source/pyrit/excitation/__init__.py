# coding=utf-8
"""
=====================================
Excitations (:mod:`pyrit.excitation`)
=====================================

.. currentmodule:: pyrit.excitation

The module Excitations handles all kinds of excitations for a problem. In the following, it is distinguished between
classes that are typically used by an user, such as those that can be applied to physical groups, i.e. that can be
seen as an actual excitation, or the container class (:py:class:`~pyrit.excitation.Excitations`) and other (abstract)
classes that are used internally. Typical excitations are charges, charge densities, currents or current densities. In
addition to that, the field-circuit coupling models solid and stranded conductor are also implemented in this package.
Like in the previous package, there is also a container class that manages all excitations of a problem.

Excitations and container class
-------------------------------

.. autosummary::
    :toctree: autosummary/

    Excitations
    ChargeDensity
    CurrentDensity
    SolidConductor
    StrandedConductor
    DistributedSolidConductor


Abstract classes
----------------

.. autosummary::
    :toctree: autosummary/

    Exci
    FieldCircuitCoupling

"""

from .Exci import Exci
from .Excitations import Excitations
from .FieldCircuitCoupling import FieldCircuitCoupling
from .SolidConductor import SolidConductor
from .StrandedConductor import StrandedConductor
from .DistributedSolidConductor import DistributedSolidConductor

from .ChargeDensity import ChargeDensity
from .CurrentDensity import CurrentDensity


__all__ = ['Exci', 'Excitations', 'FieldCircuitCoupling', 'StrandedConductor', 'SolidConductor',
           'DistributedSolidConductor', 'ChargeDensity', 'CurrentDensity']
