# coding=utf-8
"""Current flow solutions in two-dimensional Cartesian coordinates

.. sectionauthor:: Bergfried, Bundschuh
"""

from pyrit import get_logger

from pyrit.solution.Solution import StaticSolution, HarmonicSolution, TransientSolution

logger = get_logger(__name__)

__all__ = ['CurrentFlowSolutionCartStatic', 'CurrentFlowSolutionCartHarmonic', 'CurrentFlowSolutionCartTransient']


class CurrentFlowSolutionCartStatic(StaticSolution):
    """The solution of a two-dimensional stationary current problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemCartStatic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class CurrentFlowSolutionCartHarmonic(HarmonicSolution):
    """The solution of a two-dimensional, harmonic electroquasistatic problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemCartHarmonic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class CurrentFlowSolutionCartTransient(TransientSolution):
    """The solution of a two-dimensional, transient electroquasistatic problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemCartTransient`.
    """

    related_static_solution = CurrentFlowSolutionCartStatic

    def consistency_check(self):
        raise NotImplementedError()
