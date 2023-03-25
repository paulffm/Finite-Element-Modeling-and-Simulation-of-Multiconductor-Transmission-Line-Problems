# coding=utf-8
"""
================================
Solution (:mod:`pyrit.solution`)
================================

.. currentmodule:: pyrit.solution

This module is responsible for managing a problem's solution.

This module contains a variety of solution classes. They are categorized in the same was as the problem classes are.
There are the categories :ref:`label-electric-solutions`, :ref:`label-magnetic-solutions`,
:ref:`label-thermal-solutions` and :ref:`label-current-flow-solutions`.

Furthermore, there are :ref:`label-abstract-solution-classes`, :ref:`label-field-plotter` and
:ref:`label-miscellaneous`.


.. _label-electric-solutions:

Electric solutions
------------------

.. autosummary::
    :toctree: autosummary/

    ElectricSolutionCartStatic
    ElectricSolutionAxiStatic

.. _label-magnetic-solutions:

Magnetic solutions
------------------

.. autosummary::
    :toctree: autosummary/

    MagneticSolutionCartStatic
    MagneticSolutionCartHarmonic
    MagneticSolutionCartTransient
    MagneticSolutionAxiStatic
    MagneticSolutionAxiHarmonic
    MagneticSolutionAxiTransient

.. _label-thermal-solutions:

Thermal solutions
-----------------

.. autosummary::
    :toctree: autosummary/

    ThermalSolutionCartStatic
    ThermalSolutionCartHarmonic
    ThermalSolutionCartTransient
    ThermalSolutionAxiStatic
    ThermalSolutionAxiHarmonic
    ThermalSolutionAxiTransient

.. _label-current-flow-solutions:

Current flow solutions
----------------------

.. autosummary::
    :toctree: autosummary/

    CurrentFlowSolutionCartStatic
    CurrentFlowSolutionCartHarmonic
    CurrentFlowSolutionCartTransient
    CurrentFlowSolutionAxiStatic
    CurrentFlowSolutionAxiHarmonic
    CurrentFlowSolutionAxiTransient

.. _label-abstract-solution-classes:

Abstract solution classes
--------------------------------------------

.. autosummary::
    :toctree: autosummary/

    Solution
    StaticSolution
    HarmonicSolution
    TransientSolution

.. _label-field-plotter:

Field plotter
-------------
Classes for plotting static, harmonic and time dependent fields on triangular meshes.

.. autosummary::
    :toctree: autosummary/

    FieldPlotter
    StaticFieldPlotter
    HarmonicFieldPlotter
    TransientFieldPlotter

.. _label-miscellaneous:

Miscellaneous
-------------

.. autosummary::
    :toctree: autosummary/

    VectorConverter
"""

from .Solution import Solution, StaticSolution, HarmonicSolution, TransientSolution, get_field, set_field, get_matrix, \
    set_matrix

from .FieldPlotter import FieldPlotter, StaticFieldPlotter, HarmonicFieldPlotter, TransientFieldPlotter

from .ElectricSolutionCart import ElectricSolutionCartStatic
from .ElectricSolutionAxi import ElectricSolutionAxiStatic

from .MagneticSolutionCart import MagneticSolutionCartStatic, MagneticSolutionCartHarmonic, \
    MagneticSolutionCartTransient
from .MagneticSolutionAxi import MagneticSolutionAxiBase, MagneticSolutionAxiStatic, MagneticSolutionAxiHarmonic, \
    MagneticSolutionAxiTransient, VectorConverter

from .ThermalSolutionCart import ThermalSolutionCartStatic, ThermalSolutionCartHarmonic, ThermalSolutionCartTransient
from .ThermalSolutionAxi import ThermalSolutionAxiStatic, ThermalSolutionAxiHarmonic, ThermalSolutionAxiTransient

from .CurrentFlowSolutionCart import CurrentFlowSolutionCartStatic, CurrentFlowSolutionCartHarmonic, \
    CurrentFlowSolutionCartTransient
from .CurrentFlowSolutionAxi import CurrentFlowSolutionAxiStatic, CurrentFlowSolutionAxiHarmonic, \
    CurrentFlowSolutionAxiTransient

__all__ = ['Solution', 'StaticSolution', 'HarmonicSolution', 'TransientSolution', 'get_field', 'set_field',
           'get_matrix', 'set_matrix',
           'FieldPlotter', 'StaticFieldPlotter', 'HarmonicFieldPlotter', 'TransientFieldPlotter',
           'ElectricSolutionCartStatic', 'ElectricSolutionAxiStatic',
           'MagneticSolutionCartStatic', 'MagneticSolutionCartHarmonic', 'MagneticSolutionCartTransient',
           'MagneticSolutionAxiBase', 'MagneticSolutionAxiStatic', 'MagneticSolutionAxiHarmonic',
           'ThermalSolutionCartStatic', 'ThermalSolutionCartHarmonic', 'ThermalSolutionCartTransient',
           'ThermalSolutionAxiStatic', 'ThermalSolutionAxiHarmonic',
           'ThermalSolutionAxiTransient',
           'MagneticSolutionAxiTransient', 'VectorConverter',
           'CurrentFlowSolutionCartStatic', 'CurrentFlowSolutionCartHarmonic',
           'CurrentFlowSolutionCartTransient',
           'CurrentFlowSolutionAxiStatic', 'CurrentFlowSolutionAxiHarmonic',
           'CurrentFlowSolutionAxiTransient']
