# coding=utf-8
"""Two-dimensional, electrostatic problem.

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING, Tuple, Type
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy import sparse

from pyrit import get_logger
from pyrit import region, material, mesh, shapefunction


from pyrit.solution.Solution import StaticSolution
from pyrit.solution.FieldPlotter import StaticFieldPlotter

if TYPE_CHECKING:
    from pyrit import excitation
    from pyrit.shapefunction import TriCartesianNodalShapeFunction

logger = get_logger(__name__)

__all__ = ['ElectricSolutionCartStatic']


class ElectricSolutionCartStatic(StaticSolution):
    """The solution of an electrostatic problem in 2D Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ElectricProblemCartStatic`.
    """

    solution_identifier: str = 'Two-dimensional, electrostatic solution'

    ignore_for_saving = ['_divgrad_matrix']

    def __init__(self, description: str, potential: np.ndarray, tri_mesh: 'mesh.TriMesh',
                 tri_cart_nodal_shape_function: 'shapefunction.TriCartesianNodalShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: 'excitation.Excitations' = None):
        """Base class for the solution of static and harmonic two-dimensional electric problems.

        Parameters
        ----------
        description : str
            Description of the solution.
        potential : np.ndarray
            An array that represents the potential. Needs to match the number of nodes in the mesh.
        tri_mesh : mesh.TriMesh
            The mesh object.
        tri_cart_nodal_shape_function : shapefunction.TriCartesianNodalShapeFunction
            The shape function object.
        regions : region.Regions
            The regions object.
        materials : material.Materials
            The materials object.
        excitations : excitation.Excitations, optional
            The excitation object.
        """
        super().__init__(description, potential, None, None, regions, materials, excitations)

        self.mesh: mesh.TriMesh = tri_mesh
        self.shape_function: shapefunction.TriCartesianNodalShapeFunction = tri_cart_nodal_shape_function

        # The field plotter
        self.plotter = StaticFieldPlotter(tri_mesh)

        # Instantiate attributes for properties
        self._e_field = None
        self._d_field = None
        self._divgrad_matrix = None

        self.consistency_check()

    def consistency_check(self):
        super().consistency_check()

        logger.info("Starting with solution specitic consistency check on solution '%s'", str(self))

        # Type checking of types of mesh and shape function
        if not isinstance(self.mesh, mesh.TriMesh):
            logger.warning("The mesh is not of the expected type. Is %s but should be TriMesh", str(self.mesh))
        if not isinstance(self.shape_function, shapefunction.TriCartesianNodalShapeFunction):
            logger.warning("The shape function is not of the expected type. Is %s but should be "
                           "TriCartesianNodalShapeFunction", self.shape_function)

        # Check if the potential is of correct format
        self._check_potential()

    def _check_potential(self):
        """Check if the potential is of correct format"""
        if not isinstance(self.potential, np.ndarray):
            logger.warning("The potential is not a numpy array. This may lead to errors later.")
        if not self.potential.size == self.mesh.num_node:
            raise ValueError("The passed potential has a wrong number entries.")
        if not self.potential.ndim == 1:
            logger.warning("The number of dimensions of the potential is wrong. I will run a flatten on the array.")
            self.potential = self.potential.flatten()

    @property
    def divgrad_matrix(self) -> sparse.coo_matrix:
        """The divgrad matrix."""
        if self._divgrad_matrix is None:
            self._divgrad_matrix = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                        material.Permittivity)
        return self._divgrad_matrix

    @divgrad_matrix.setter
    def divgrad_matrix(self, divgrad_matrix: sparse.coo_matrix):
        self._divgrad_matrix = divgrad_matrix

    @property
    def potential(self) -> np.ndarray:
        """The potential"""
        return self._solution

    @potential.setter
    def potential(self, potential: np.ndarray):
        self._solution = potential
        self._check_potential()
        self._e_field = None
        self._d_field = None

    @property
    def e_field(self) -> np.ndarray:
        """The electric field strength."""
        if self._e_field is None:
            self._e_field = -1 * self.shape_function.gradient(self.potential)
        return self._e_field

    @e_field.setter
    def e_field(self, e_field: np.ndarray):
        self._e_field = e_field

    @property
    def d_field(self) -> np.ndarray:
        """The displacement current."""
        if self._d_field is None:
            permittivity_per_element = self.mesh.value_per_element(self.regions, self.materials, material.Permittivity)
            self._d_field = np.tile(np.atleast_2d(permittivity_per_element).T, 2) * self.e_field
        return self._d_field

    @d_field.setter
    def d_field(self, d_field: np.ndarray):
        self._d_field = d_field

    @property
    def energy_density(self):
        """The energy density"""
        return 0.5 * np.sum(self.e_field * self.d_field, axis=1)

    @property
    def energy(self):
        """The energy in the whole domain."""
        return 0.5 * self.potential @ self.divgrad_matrix @ self.potential

    def energy_in_regions(self, *region_ids: int) -> float:
        """Compute the energy in certain regions.

        Parameters
        ----------
        region_ids : int
            Region IDs of the regions to compute the energy in.

        Returns
        -------
        energy : float
            The energy

        Raises
        ------
        ValueError : When the region of one ID is not two-dimensional.
        """
        indices_per_region = []  # List of indices of the elements in the regions
        for region_id in region_ids:
            if not self.regions.get_regi(region_id).dim == 2:
                raise ValueError(f"The region '{region_id}' does not belong to a region of the expected dimension.")
            indices_per_region.append(np.where(self.mesh.elem2regi == region_id))

        indices = np.concatenate(indices_per_region)  # List of all indices in the regions

        energy = self.shape_function.length * np.sum(self.energy_density[indices] * self.mesh.elem_area[indices])

        return energy

    def plot_potential(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the potential.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`StaticFieldPlotter.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Potential in V'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_scalar_field(self.potential, **kwargs)

    def plot_equilines(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot equilines for the potential.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`StaticFieldPlotter.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Equipotential lines'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_equilines(self.potential, **kwargs)

    def plot_e_field(self, plot_type: str, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the electric field strength.

        Parameters
        ----------
        plot_type : str
            The plot type of the plot. See :py:meth:`StaticFieldPlotter.plot_scalar_field`.
        kwargs :
            Passed to :py:meth:`StaticFieldPlotter.plot_vector_field`.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Electric field strength $\vec{E}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.e_field, plot_type, **kwargs)

    def plot_d_field(self, plot_type: str, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the displacement current.

        Parameters
        ----------
        plot_type : str
            The plot type of the plot. See :py:meth:`StaticFieldPlotter.plot_scalar_field`.
        kwargs :
            Passed to :py:meth:`StaticFieldPlotter.plot_vector_field`.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Displacement current $\vec{D}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.d_field, plot_type, **kwargs)

    def plot_energy_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the energy density.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`StaticFieldPlotter.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Energy density'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.energy_density, **default_kwargs)
