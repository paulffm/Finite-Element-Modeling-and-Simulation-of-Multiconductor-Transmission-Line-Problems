# coding=utf-8
"""Thermal problems in two-dimensional Cartesian coordinates

.. sectionauthor:: Bundschuh
"""
from typing import Dict, Union, Tuple, Callable, NoReturn

import numpy as np

from pyrit.solution import ThermalSolutionCartStatic, ThermalSolutionCartHarmonic, ThermalSolutionCartTransient
from pyrit.problem.Problem import StaticProblem, HarmonicProblem, TransientProblem, Solution_Monitor, Monitor


class ThermalProblemCartStatic(StaticProblem):
    r"""A two-dimensional stationary heat conduction problem in Cartesian coordinates:

    The stationary heat conduction problem models the steady-state temperature distribution of a body with \
    constant external heat sources. The corresponding differential equation reads

    .. math::
        -\mathrm{div}(\lambda \, \mathrm{grad} (T)) = \dot q,

    where :math:`\lambda` is the thermal conductivity (see :py:class:`ThermalConductivity`), \
    :math:`\dot q` is the heat flux density and :math:`T` denotes the temperature.
    """

    problem_identifier: str = 'Stationary heat conduction problem in two-dimensional Cartesian coordinates'

    def solve(self, *args, **kwargs) -> ThermalSolutionCartStatic:
        raise NotImplementedError()


class ThermalProblemCartHarmonic(HarmonicProblem):
    r"""A two-dimensional harmonic heat conduction problem in Cartesian coordinates:

    The harmonic heat conduction problem models the steady-state temperature distribution of a body with \
    harmonic external heat sources. The corresponding differential equation reads

        .. math::
            j \omega s T -\mathrm{div}(\lambda\,\mathrm{grad}(T)) = \dot q,

    where :math:`\lambda` is the thermal conductivity (see :py:class:`ThermalConductivity`), \
    :math:`s` is the volumetric heat capacity (see :py:class:`VolumetricHeatCapacity`), \
    :math:`\dot q` is the heat flux density and :math:`T` denotes the temperature.
    """

    problem_identifier: str = 'Harmonic heat conduction problem in two-dimensional Cartesian coordinates'

    def solve(self, *args, **kwargs) -> ThermalSolutionCartHarmonic:
        raise NotImplementedError()


class ThermalProblemCartTransient(TransientProblem):
    r"""A two-dimensional transient heat conduction problem in Cartesian coordinates:

    This problem models the transient conductive heat transfer inside a body. The corresponding differential \
    equation reads

        .. math::
            \partial_t (s T) -\mathrm{div}(\lambda\,\mathrm{grad}(T)) = \dot q,

    where :math:`\lambda` is the thermal conductivity (see :py:class:`ThermalConductivity`), \
    :math:`s` is the volumetric heat capacity (see :py:class:`VolumetricHeatCapacity`), \
    :math:`\dot q` is the heat flux density and :math:`T` denotes the temperature.
    """

    problem_identifier: str = 'Transient heat conduction problem in two-dimensional Cartesian coordinates'

    def solve(self, start_value: np.ndarray, solution_monitor: Solution_Monitor = 1,
              monitors: Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]] = None,
              callback: Callable[['ThermalSolutionCartStatic'], NoReturn] = None,
              **kwargs) -> ThermalSolutionCartTransient:
        raise NotImplementedError()
