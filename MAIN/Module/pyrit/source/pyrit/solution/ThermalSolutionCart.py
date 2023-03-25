# coding=utf-8
"""Thermal solutions in two-dimensional Cartesian coordinates

.. sectionauthor:: Bundschuh
"""

from pyrit import get_logger

from pyrit.solution.Solution import StaticSolution, HarmonicSolution, TransientSolution

logger = get_logger(__name__)


class ThermalSolutionCartStatic(StaticSolution):
    """The solution of a two-dimensional, stationary heat conduction problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemCartStatic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class ThermalSolutionCartHarmonic(HarmonicSolution):
    """The solution of a two-dimensional, harmonic heat conduction problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemCartHarmonic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class ThermalSolutionCartTransient(TransientSolution):
    """The solution of a two-dimensional, transient heat conduction problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemCartTransient`.
    """

    related_static_solution = ThermalSolutionCartStatic

    def plot_monitor_solution(self, key: str, **kwargs):
        raise NotImplementedError()

    def consistency_check(self):
        raise NotImplementedError()
