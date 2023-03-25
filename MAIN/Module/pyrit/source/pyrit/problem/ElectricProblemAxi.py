# coding=utf-8
"""Electric problems in axisymmetric coordinates

.. sectionauthor:: Bundschuh
"""

from pyrit.solution import ElectricSolutionAxiStatic
from pyrit.problem.Problem import StaticProblem


class ElectricProblemAxiStatic(StaticProblem):
    r"""An electrostatic problem in axisymmetric coordinates:

    The electrostatic problem models capacitive effects, the corresponding differential equation reads

    .. math::
        -\mathrm{div}(\varepsilon \, \mathrm{grad} (\phi)) = \varrho,

    where :math:`\varepsilon` is the electric permittivity (see :py:class:`Permittivity`), :math:`\phi` is the \
    electric scalar potential and :math:`\varrho` is the charge density. A possible application is, for example, the \
    electric field caused by a point charge in space.
    """

    problem_identifier: str = 'Electrostatic problem in axisymmetric coordinates'

    def solve(self, *args, **kwargs) -> ElectricSolutionAxiStatic:
        raise NotImplementedError()
