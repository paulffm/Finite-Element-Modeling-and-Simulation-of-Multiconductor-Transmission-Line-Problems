# coding=utf-8
"""Electric solutions in axisymmetric coordinates

.. sectionauthor:: Bundschuh
"""

from pyrit import get_logger

from pyrit.solution.Solution import StaticSolution

logger = get_logger(__name__)


class ElectricSolutionAxiStatic(StaticSolution):
    """The solution of an electrostatic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ElectricProblemAxiStatic`.
    """

    def consistency_check(self):
        raise NotImplementedError()
