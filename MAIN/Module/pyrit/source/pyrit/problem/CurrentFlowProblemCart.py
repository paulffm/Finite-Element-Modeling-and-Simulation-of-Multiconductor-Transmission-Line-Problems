# coding=utf-8
"""Current flow problem in two-dimensional Cartesian coordinates

.. sectionauthor:: Bergfried, Bundschuh
"""
from typing import Dict, Union, Tuple, Callable, NoReturn

import numpy as np

from pyrit import get_logger

from pyrit.solution import CurrentFlowSolutionCartStatic, CurrentFlowSolutionCartHarmonic, \
    CurrentFlowSolutionCartTransient
from pyrit.problem.Problem import StaticProblem, HarmonicProblem, TransientProblem, Solution_Monitor, Monitor


logger = get_logger(__name__)

__all__ = ['CurrentFlowProblemCartStatic', 'CurrentFlowProblemCartHarmonic', 'CurrentFlowProblemCartTransient']


class CurrentFlowProblemCartStatic(StaticProblem):
    r"""A two-dimensional stationary current problem in Cartesian coordinates:

    The stationary current problem models resistive effects. The corresponding differential equation reads

    .. math::
        -\mathrm{div}(\sigma \, \mathrm{grad} (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`Conductivity`) and :math:`\phi` denotes the \
    electric scalar potential. A possible application is, for example, the steady state simulation of a HVDC cable \
    joint.
    """

    problem_identifier: str = 'Static, current flow problem in two-dimensional Cartesian coordinates'

    def solve(self, *args, **kwargs) -> CurrentFlowSolutionCartStatic:
        raise NotImplementedError()


class CurrentFlowProblemCartHarmonic(HarmonicProblem):
    r"""A two-dimensional harmonic electroquasistatic problem in Cartesian coordinates:

    The harmonic electroquasistatic problem models capacitive-resistive effects. The corresponding
    differential \
    equation reads

    .. math::
        -\mathrm{div}(j 2\pi f\varepsilon\, \mathrm{grad} (\phi)) -\mathrm{div}(\sigma\,\mathrm{grad}
         (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`~pyrit.material.Conductivity`),
    :math:`\varepsilon` is the electric permittivity (see :py:class:`~pyrit.material.Permittivity`), :math:`\phi`
    denotes the electric scalar potential, and :math:`f` is the frequency. The problem can, by definition,
    only handle linear materials.
    A possible application is, e.g. the steady state simulation of the insulation system in an electrical machine.
    """

    problem_identifier = 'Harmonic, electroquasistatic problem in 2D Cartesian coordinates'

    def solve(self, *args, **kwargs) -> CurrentFlowSolutionCartHarmonic:
        raise NotImplementedError()


class CurrentFlowProblemCartTransient(TransientProblem):
    r"""A two-dimensional transient electroquasistatic problem in Cartesian coordinates:

    The transient electroquasistatic problem models capacitive-resistive effects. The corresponding differential \
    equation reads

    .. math::
        -\mathrm{div}(\partial_t(\varepsilon\, \mathrm{grad} (\phi))) -\mathrm{div}(\sigma\,\mathrm{grad}
         (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`Conductivity`), :math:`\varepsilon` is the \
    electric permittivity (see :py:class:`Permittivity`) and :math:`\phi` denotes the electric scalar potential. \
    A possible application is, for example, the simulation of transient effects in HVDC cable joints or AC surge \
    arresters.
    """

    problem_identifier = 'Transient, electroquasistatic problem in 2D Cartesian coordinates'

    def solve(self, start_value: np.ndarray, solution_monitor: Solution_Monitor = 1,
              monitors: Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]] = None,
              callback: Callable[['CurrentFlowProblemCartStatic'], NoReturn] = None,
              **kwargs) -> CurrentFlowSolutionCartTransient:
        raise NotImplementedError()
