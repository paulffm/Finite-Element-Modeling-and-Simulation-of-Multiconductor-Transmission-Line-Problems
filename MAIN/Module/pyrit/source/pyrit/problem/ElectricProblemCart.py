# coding=utf-8
"""Two-dimensional, electrostatic problem.

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING

from pyrit import get_logger
from pyrit import mesh, shapefunction, region, material, bdrycond, excitation

from pyrit.solution import ElectricSolutionCartStatic
from pyrit.problem.Problem import StaticProblem

if TYPE_CHECKING:
    from pyrit.mesh import TriMesh
    from pyrit.shapefunction import TriCartesianNodalShapeFunction

logger = get_logger(__name__)

__all__ = ['ElectricProblemCartStatic']


class ElectricProblemCartStatic(StaticProblem):
    r"""A two-dimensional, electrostatic problem in Cartesian coordinates:

    The electrostatic problem models capacitive effects, the corresponding differential equation reads

    .. math::
        -\mathrm{div}(\varepsilon \, \mathrm{grad} (\phi)) = \varrho\,,

    where :math:`\varepsilon` is the electric permittivity (see :py:class:`~pyrit.material.Permittivity`), :math:`\phi`
    is the electric scalar potential and :math:`\varrho` is the charge density. A possible application is, for example,
    the electric field caused by a point charge in space.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.ElectricSolutionCartStatic`.
    """

    problem_identifier = 'Electrostatic problem in 2D Cartesian coordinates'

    def __init__(self, description: str, tri_mesh: mesh.TriMesh,
                 tri_cart_nodal_shape_function: shapefunction.TriCartesianNodalShapeFunction,
                 regions: region.Regions, materials: material.Materials, boundary_conditions: bdrycond.BdryCond,
                 excitations: excitation.Excitations = None):
        """An electrostatic problem in 2D Cartesian coordinates.

        Parameters
        ----------
        description : str
            A description of the problem
        tri_mesh : TriMesh
            A mesh object.
        tri_cart_nodal_shape_function : TriCartesianNodalShapeFunction
            A shape function object.
        regions : Regions
            A regions object.
        materials : Materials
            A materials object.
        boundary_conditions : BdryCond
            A boundary conditions object.
        excitations : Excitations
            An excitations object.
        """
        super().__init__(description, None, None, regions, materials, boundary_conditions, excitations)

        self.mesh: mesh.TriMesh = tri_mesh
        self.shape_function: shapefunction.TriCartesianNodalShapeFunction = tri_cart_nodal_shape_function

        self.consistency_check()

    def consistency_check(self):
        super().consistency_check()

        logger.info("Starting with problem specific consistency check on problem '%s'.", str(self))
        # Type checking of types of mesh and shape function
        if not isinstance(self.mesh, mesh.TriMesh):
            logger.warning("The mesh is not of the expected type. Is %s but should be TriMesh", str(self.mesh))
        if not isinstance(self.shape_function, shapefunction.TriCartesianNodalShapeFunction):
            logger.warning("The shape function is not of the expected type. Is %s but should be "
                           "TriCartesianNodalShapeFunction", self.shape_function)

        # Check if all regions are in the mesh
        self._check_all_regions_in_trimesh()
        logger.info("Done with problem specific consistency check on problem '%s'.", self.description)

    def solve(self, *args, **kwargs) -> ElectricSolutionCartStatic:
        """Solve the problem.

        Parameters
        ----------
        args :
        kwargs :
            Passed to :py:func:`solve_linear_system`.

        Returns
        -------
        solution : ElectricSolutionCartStatic
        """
        divgrad = self.shape_function.divgrad_operator(self.regions, self.materials, material.Permittivity)
        if self.excitations is None:
            q = self.shape_function.load_vector(0)
        else:
            q = self.shape_function.load_vector(self.regions, self.excitations)

        matrix_shrink, q_shrink, _, _, support_data = self.shape_function.shrink(divgrad, q, self, 1)
        u_shrink, _ = type(self).solve_linear_system(matrix_shrink.tocsr(), q_shrink.tocsr(), **kwargs)
        potential = self.shape_function.inflate(u_shrink, self, support_data)

        solution = ElectricSolutionCartStatic(f'solution to self {self.description}', potential,
                                              self.mesh, self.shape_function, self.regions,
                                              self.materials)
        solution.divgrad_matrix = divgrad

        return solution
