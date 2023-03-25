# coding=utf-8
"""
===========================================
Boundary conditions (:mod:`pyrit.bdrycond`)
===========================================

.. currentmodule:: pyrit.bdrycond

This modules defines the supported boundary conditions. Besides the standard Dirichlet and Neumann boundary
conditions, also Robin, binary and floating boundary conditions are defined. There is a container class that manages
all boundary conditions of a problem and can apply them to the system of equations of the corresponding problem.

Boundary conditions and container class
---------------------------------------

.. autosummary::
    :toctree: autosummary/

    BdryCond
    BCAntiPeriodic
    BCBinary
    BCDirichlet
    BCDummy
    BCFloating
    BCNeumann
    BCPeriodic
    BCRobin

Abstract classes
----------------

.. autosummary::
    :toctree: autosummary/

    BC
    BCUnary
"""

from .BC import BC
from .BCUnary import BCUnary
from .BCBinary import BCBinary
from .BCDirichlet import BCDirichlet
from .BCDummy import BCDummy
from .BCFloating import BCFloating
from .BCNeumann import BCNeumann
from .BCPeriodic import BCPeriodic
from .BCRobin import BCRobin
from .BCAntiPeriodic import BCAntiPeriodic
from .BdryCond import BdryCond


__all__ = ['BC', 'BCAntiPeriodic', 'BCBinary', 'BCDirichlet', 'BCDummy', 'BCFloating', 'BCNeumann', 'BCPeriodic',
           'BCRobin', 'BCUnary', 'BdryCond']
