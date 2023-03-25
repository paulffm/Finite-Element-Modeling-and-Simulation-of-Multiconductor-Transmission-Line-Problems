# coding=utf-8
"""Two-dimensional magnetic problem.

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING, Dict, Union, Tuple, Callable, NoReturn
import numpy as np

from pyrit import get_logger
from pyrit import mesh, shapefunction, region, material, bdrycond, excitation

from pyrit.solution import MagneticSolutionCartStatic, MagneticSolutionCartHarmonic, MagneticSolutionCartTransient
from pyrit.problem.Problem import StaticProblem, HarmonicProblem, TransientProblem, Solution_Monitor, Monitor

if TYPE_CHECKING:
    from pyrit.mesh import TriMesh

logger = get_logger(__name__)

__all__ = ['MagneticProblemCartStatic', 'MagneticProblemCartHarmonic', 'MagneticProblemCartTransient']


class MagneticProblemCartStatic(StaticProblem):
    r"""A static, two-dimensional magnetic problem in Cartesian coordinates.

    The corresponding differential equations reads:

    .. math::

        \mathrm{rot}(\nu\,\mathrm{rot}\vec{A}) = \vec{J}_{\mathrm{s}}\,,

    where :math:`\nu` is the reluctivity (see :py:class:`~pyrit.material.Reluctivity`), :math:`\vec{A}` is the magnetic
    vector potential and :math:`\vec{J}_\mathrm{s}` is the source current density. This problem can be used, whenever
    the static magnetic field that is caused by a current shall be computed. A very simple example is therefore the
    magnetic field in a coaxial cable caused by a constant current.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.MagneticSolutionCartStatic`.
    """

    problem_identifier = 'Magnetic 2D problem'

    def __init__(self, description: str, tri_cart_edge_shape_function: shapefunction.TriCartesianEdgeShapeFunction,
                 tri_mesh: mesh.TriMesh, regions: region.Regions, materials: material.Materials,
                 boundary_conditions: bdrycond.BdryCond, excitations: excitation.Excitations = None):
        """Constructor

        Parameters
        ----------
        description : str
            A description of the problem
        tri_mesh : TriMesh
            A mesh object. See :py:class:`pyrit.mesh.TriMesh`
        tri_cart_edge_shape_function : TriCartesianEdgeShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriCartesianEdgeShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`
        boundary_conditions : BdryCond
            A boundary conditions object. See :py:mod:`pyrit.bdrycond`
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`
        """
        super().__init__(description, None, None, regions, materials, boundary_conditions, excitations)

        self.mesh: mesh.TriMesh = tri_mesh
        self.shape_function: shapefunction.TriCartesianEdgeShapeFunction = tri_cart_edge_shape_function

        self.consistency_check()

    def consistency_check(self):
        super().consistency_check()

        logger.info("Starting with problem specific consistency check on problem '%s'.", self.description)
        # Type checking of types of mesh and shape function
        if not isinstance(self.mesh, mesh.TriMesh):
            logger.warning("The mesh is not of the expected type. Is %s but should be TriMesh", str(self.mesh))
        if not isinstance(self.shape_function, shapefunction.TriCartesianEdgeShapeFunction):
            logger.warning("The shape function is not of the expected type. Is %s but should be "
                           "TriCartesianNodalShapeFunction", self.shape_function)

        # Check if all regions are in the mesh
        self._check_all_regions_in_trimesh()

        logger.info("Done with problem specific consistency check on problem '%s'.", self.description)

    @property
    def length(self):
        """The length of the problem. Taken from the shape function object."""
        return self.shape_function.length

    def solve(self, *args, **kwargs) -> MagneticSolutionCartStatic:
        curlcurl = self.shape_function.curlcurl_operator(self.regions, self.materials, material.Reluctivity)
        load = self.shape_function.load_vector(self.regions, self.excitations)

        matrix_shrink, rhs_shrink, _, _, support_data = self.shape_function.shrink(curlcurl, load, self, 1)
        a_shrink, _ = type(self).solve_linear_system(matrix_shrink.tocsr(), rhs_shrink.tocsr(), **kwargs)
        vector_potential = self.shape_function.inflate(a_shrink, self, support_data)

        solution = MagneticSolutionCartStatic(f'Solution to problem \'{self.description}\'', vector_potential,
                                              self.mesh, self.shape_function, self.regions,
                                              self.materials, self.excitations)
        solution.curlcurl_matrix = curlcurl

        return solution


class MagneticProblemCartHarmonic(HarmonicProblem):
    r"""A harmonic, magnetic problem in two-dimensional Cartesian coordinates.

    This problem models magnetic fields caused by currents and takes eddy current effects into account.
    The corresponding differential equation is:

    .. math::

        \mathrm{rot}(\nu\,\mathrm{rot}\vec{A}) + j\omega\sigma\vec{A} = \vec{J}_{\mathrm{s}}\,,

    where :math:`\nu` is the Reluctivity (see :py:class:`~pyrit.material.Reluctivity`), :math:`\sigma` is the
    conductivity (see :py:class:`~pyrit.material.Conductivity`), :math:`\vec{A}` is the magnetic vector potential,
    :math:`\vec{J}_\mathrm{s}` is the source current density, and :math:`\omega=2\pi f` is the angular frequency, with
    the frequency :math:`f`.
    Both, the magnetic vector potential and the source current density are complex. The problem can, by definition,
    only handle linear materials.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.MagneticSolutionCartHarmonic`.
    """

    problem_identifier: str = "Harmonic, magnetic problem in two-dimensional Cartesian coordinates"

    def solve(self, *args, **kwargs) -> MagneticSolutionCartHarmonic:
        raise NotImplementedError()


class MagneticProblemCartTransient(TransientProblem):
    r"""A Transient, magnetic problem in two-dimensional Cartesian coordinates.

    This problem models magnetic fields caused by currents and takes eddy current effects into account.
    The corresponding differential equation is

    .. math::

        \mathrm{rot}(\nu\,\mathrm{rot}\vec{A}) + \sigma\frac{\partial\vec{A}}{\partial t} = \vec{J}_{\mathrm{s}}\,,

    where :math:`\nu` is the Reluctivity (see :py:class:`~pyrit.material.Reluctivity`), :math:`\sigma` is the
    conductivity (see :py:class:`~pyrit.material.Conductivity`), :math:`\vec{A}` is the magnetic vector potential,
    :math:`\vec{J}_\mathrm{s}` is the source current density, and :math:`t` is the time.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.MagneticSolutionCartTransient`.
    """

    problem_identifier: str = "Transient, magnetic problem in two-dimensional Cartesian coordinates"

    def solve(self, start_value: np.ndarray, solution_monitor: Solution_Monitor = 1,
              monitors: Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]] = None,
              callback: Callable[['MagneticProblemCartStatic'], NoReturn] = None,
              **kwargs) -> MagneticSolutionCartTransient:
        raise NotImplementedError()
