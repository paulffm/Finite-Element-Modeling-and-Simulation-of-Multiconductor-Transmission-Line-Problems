# coding=utf-8
"""
==============================
Problem (:mod:`pyrit.problem`)
==============================

.. currentmodule:: pyrit.problem

This module is responsible for managing problems.

Predefined Problem/Solution Types
---------------------------------

**Pyrit** offers classes for electromagnetic and heat transfer problems of any time dependence. This allows to solve
field problems based on a template-based workflow. Here, each supported problem type is listed and
associated with the partial differential equation that governs the problem.

This package collects classes that represent standard problem formulations. They provide a framework for organizing
all the data needed for defining a problem. Furthermore, they
implement a ``solve`` method that provides a convenient way to solve a problem and returns a solution object (see the
:py:mod:`~pyrit.solution`). Per problem type, i.e. the underlying differential equation, there are separate
classes for static, harmonic and transient problems.

Pyrit supports *magnetic, electric, current flow, and thermal modules* of problem solution counterparts. Note
that the modules :py:mod:`ElectricProblemCart` and :py:mod:`ElectricProblemAxi` solve the Poisson equation and,
therefore, only provide static classes. *Electroquasistatic* scenarios can be simulated with modules
:py:mod:`CurrentFlowProblemCart` and :py:mod:`CurrentFlowProblemAxi`.

Each problem can be computed for *static*,
*harmonic*, and *transient* scenarios.

For the underlying FE discretization ansatz, it is important to distinguish
*Cartesian* or *axisymmetric* cases. Therefore, for each case, i.e. *Cartesian* case (:math:`u(x,y)`), denoted as
:py:mod:`ProblemCart` and *axisymmetric* case (:math:`u(r,z)`) denoted as :py:mod:`ProblemAxi`, a respective problem
and solution pair exists.

The problem classes can be categorized into :ref:`label-electric-problems`, :ref:`label-magnetic-problems`,
:ref:`label-thermal-problems` and :ref:`label-current-flow-problems`.
Furthermore there are :ref:`label-general-problem-classes`.

.. _label-electric-problems:

Electric problems
-----------------

.. autosummary::
    :toctree: autosummary/

    ElectricProblemCartStatic
    ElectricProblemAxiStatic

.. _label-magnetic-problems:

Magnetic problems
-----------------

.. autosummary::
    :toctree: autosummary/

    MagneticProblemCartStatic
    MagneticProblemCartHarmonic
    MagneticProblemCartTransient
    MagneticProblemAxiStatic
    MagneticProblemAxiHarmonic
    MagneticProblemAxiTransient

.. _label-thermal-problems:

Thermal problems
----------------

.. autosummary::
    :toctree: autosummary/

    ThermalProblemCartStatic
    ThermalProblemCartHarmonic
    ThermalProblemCartTransient
    ThermalProblemAxiStatic
    ThermalProblemAxiHarmonic
    ThermalProblemAxiTransient
    ThermalSolverInfoAxiStatic
    ThermalSolverInfoAxiTransient

.. _label-current-flow-problems:

Current flow problems
----------------------

.. autosummary::
    :toctree: autosummary/

    CurrentFlowSolverInfoAxiStatic
    CurrentFlowSolverInfoAxiTransient
    CurrentFlowProblemCartStatic
    CurrentFlowProblemCartHarmonic
    CurrentFlowProblemCartTransient
    CurrentFlowProblemAxiStatic
    CurrentFlowProblemAxiHarmonic
    CurrentFlowProblemAxiTransient

.. _label-general-problem-classes:

General problem classes
-----------------------

.. autosummary::
    :toctree: autosummary/

    Problem
    StaticProblem
    HarmonicProblem
    TransientProblem
"""

from .Problem import Problem, StaticProblem, HarmonicProblem, TransientProblem

from .CurrentFlowProblemAxi import CurrentFlowProblemAxiStatic, CurrentFlowProblemAxiHarmonic, \
    CurrentFlowProblemAxiTransient, CurrentFlowSolverInfoAxiStatic, CurrentFlowSolverInfoAxiTransient
from .CurrentFlowProblemCart import CurrentFlowProblemCartStatic, CurrentFlowProblemCartHarmonic, \
    CurrentFlowProblemCartTransient

from .ElectricProblemAxi import ElectricProblemAxiStatic
from .ElectricProblemCart import ElectricProblemCartStatic

from .MagneticProblemAxi import MagneticProblemAxiStatic, MagneticProblemAxiHarmonic, MagneticProblemAxiTransient
from .MagneticProblemCart import MagneticProblemCartStatic, MagneticProblemCartHarmonic, MagneticProblemCartTransient

from .ThermalProblemAxi import ThermalProblemAxiStatic, ThermalProblemAxiHarmonic, ThermalProblemAxiTransient, \
    ThermalSolverInfoAxiStatic, ThermalSolverInfoAxiTransient
from .ThermalProblemCart import ThermalProblemCartStatic, ThermalProblemCartHarmonic, ThermalProblemCartTransient

__all__ = ['Problem', 'StaticProblem', 'HarmonicProblem', 'TransientProblem',
           'CurrentFlowProblemAxiStatic', 'CurrentFlowProblemAxiHarmonic', 'CurrentFlowProblemAxiTransient',
           'CurrentFlowSolverInfoAxiStatic', 'CurrentFlowSolverInfoAxiTransient',
           'CurrentFlowProblemCartStatic', 'CurrentFlowProblemCartHarmonic', 'CurrentFlowProblemCartTransient',
           'ElectricProblemAxiStatic',
           'ElectricProblemCartStatic',
           'MagneticProblemAxiStatic', 'MagneticProblemAxiHarmonic', 'MagneticProblemAxiTransient',
           'MagneticProblemCartStatic', 'MagneticProblemCartHarmonic', 'MagneticProblemCartTransient',
           'ThermalProblemAxiStatic', 'ThermalProblemAxiHarmonic', 'ThermalProblemAxiTransient',
           'ThermalSolverInfoAxiStatic', 'ThermalSolverInfoAxiStatic',
           'ThermalProblemCartStatic', 'ThermalProblemCartHarmonic', 'ThermalProblemCartTransient']
