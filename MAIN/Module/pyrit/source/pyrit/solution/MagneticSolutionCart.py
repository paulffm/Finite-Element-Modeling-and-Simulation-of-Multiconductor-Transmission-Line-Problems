# coding=utf-8
"""Magnetic solution Cartesian

.. sectionauthor:: Bundschuh
"""

from typing import TYPE_CHECKING, Tuple, Type

import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy.sparse import linalg as splinalg
from scipy import sparse

from pyrit import get_logger
from pyrit import region, material
from pyrit import excitation

from pyrit.solution.Solution import StaticSolution, HarmonicSolution, TransientSolution
from pyrit.solution.FieldPlotter import StaticFieldPlotter

if TYPE_CHECKING:
    from pyrit import shapefunction, mesh
    from pyrit.shapefunction import TriCartesianEdgeShapeFunction

logger = get_logger(__name__)


# noinspection PyAttributeOutsideInit
class MagneticSolutionCartStatic(StaticSolution):
    r"""Solution of a static, two-dimensional magnetic problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemCartStatic`.
    """

    solution_identifier = 'Solution of a static, tow-dimensional magnetostatic problem in Cartesian coordinates'

    #: Protected attribute names that are used for properties
    _property_attributes = ['_b_field', '_h_field', '_curlcurl_matrix', '_mass_one', '_load_vector']
    ignore_for_saving = _property_attributes.copy()

    def __init__(self, description: str, vector_potential: np.ndarray, tri_mesh: 'mesh.TriMesh',
                 tri_cart_edge_shape_function: 'shapefunction.TriCartesianEdgeShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations):
        """Constructor

        Parameters
        ----------
        description : str
            Description of the solution.
        vector_potential : np.ndarray
            An array that represents the potential.
        tri_mesh : mesh.TriMesh
            The mesh object.
        tri_cart_edge_shape_function : shapefunction.TriCartesianEdgeShapeFunction
            The shape function object.
        regions : region.Regions
            The regions object.
        materials : material.Materials
            The materials object.
        excitations : excitations.Excitations
            The excitations object.
        """
        super().__init__(description, vector_potential, None, None, regions, materials, excitations)

        self.mesh: 'mesh.TriMesh' = tri_mesh
        self.shape_function: 'shapefunction.TriCartesianEdgeShapeFunction' = tri_cart_edge_shape_function

        self.plotter = StaticFieldPlotter(tri_mesh)

        self.consistency_check()

    def consistency_check(self):
        # Check if the vector_potential is of the right format
        if not self.vector_potential.ndim == 1:
            logger.warning("The vector-potential is not of expected format. It is not 1 dimensional")
            if self.vector_potential.ndim == 2:
                if self.vector_potential.size == self.mesh.num_node:
                    logger.info("I will flatten the vector potential")
                    self.vector_potential = self.vector_potential.flatten()

        if not issubclass(self.vector_potential.dtype.type, (np.integer, np.floating)):
            logger.warning("The vector potential is of unexpected type.")

    @property
    def length(self) -> float:
        """The length of the domain."""
        return self.shape_function.length

    @property
    def vector_potential(self) -> np.ndarray:
        r"""The vector potential. These are the coefficients for the basis functions.

        Note that this vector does not represent the values of :math:`\vec{A}`. See the property
        vector_potential_field_values.
        """
        return self._solution

    @vector_potential.setter
    def vector_potential(self, vector_potential: np.ndarray):
        self._solution = vector_potential

    @property
    def vector_potential_field_values(self):
        """The field values of the vector potential."""
        return self.vector_potential / self.length

    @property
    def b_field(self) -> np.ndarray:
        """Magnetic flux density. Array with two values (x- and y-component) per mesh element."""
        if self._b_field is None:
            self._b_field = self.shape_function.curl(self.vector_potential)
        return self._b_field

    @b_field.setter
    def b_field(self, b_field: np.ndarray):
        self._b_field = b_field

    @property
    def h_field(self) -> np.ndarray:
        """Magnetic field strength. Array with two values (x- and y-component) per mesh element."""
        if self._h_field is None:
            reluctivity_per_element = self.mesh.value_per_element(self.regions, self.materials, material.Reluctivity)
            self._h_field = np.tile(np.atleast_2d(reluctivity_per_element).T, 2) * self.b_field
        return self._h_field

    @h_field.setter
    def h_field(self, h_field: np.ndarray):
        self._h_field = h_field

    @property
    def curlcurl_matrix(self) -> sparse.coo_matrix:
        """The curl-curl-matrix."""
        if self._curlcurl_matrix is None:
            self._curlcurl_matrix = self.shape_function.curlcurl_operator(self.regions, self.materials,
                                                                          material.Reluctivity)
        return self._curlcurl_matrix

    @curlcurl_matrix.setter
    def curlcurl_matrix(self, curlcurl_matrix: sparse.coo_matrix):
        self._curlcurl_matrix = curlcurl_matrix

    @property
    def load_vector(self) -> np.ndarray:
        """The load vector. Array with one value per mesh node."""
        if self._load_vector is None:
            self._load_vector = self.shape_function.load_vector(self.regions, self.excitations).toarray()
        return self._load_vector

    @load_vector.setter
    def load_vector(self, load_vector: np.ndarray):
        self._load_vector = load_vector

    @property
    def mass_one(self) -> sparse.coo_matrix:
        """The identity mass matrix."""
        if self._mass_one is None:
            self._mass_one = self.shape_function.mass_matrix(1)
        return self._mass_one

    @mass_one.setter
    def mass_one(self, mass_one: sparse.coo_matrix):
        self._mass_one = mass_one

    @property
    def current_density(self) -> np.ndarray:
        """The current density."""
        return splinalg.spsolve(self.mass_one.tocsr(), self.load_vector)

    @property
    def energy(self) -> float:
        """The energy in the whole domain."""
        return 0.5 * self.vector_potential @ self.curlcurl_matrix @ self.vector_potential

    @property
    def energy_density(self) -> np.ndarray:
        """The energy density. Array with one value per mesh node."""
        return 0.5 * np.sum(self.b_field * self.h_field, axis=1)

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
            indices_per_region.append(np.where(self.mesh.elem2regi == region_id)[0])

        indices = np.concatenate(indices_per_region)  # List of all indices in the regions

        energy = self.length * np.sum(self.energy_density[indices] * self.mesh.elem_area[indices])

        return energy

    def plot_vector_potential(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the vector potential.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Vector potential'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.vector_potential_field_values, **kwargs)

    def plot_b_field(self, plot_type: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the magnetic flux density.

        Parameters
        ----------
        plot_type : str
            The plot_type of the plot. See :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_vector_field`
        kwargs :
            Passed to the plot method.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic flux density $\vec{B}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.b_field, plot_type, **kwargs)

    def plot_h_field(self, plot_type: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the magnetic field strength.

        Parameters
        ----------
        plot_type : str
            The plot_type of the plot. See :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_vector_field`
        kwargs :
            Passed to the plot method.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic field strength $\vec{H}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.h_field, plot_type, **kwargs)

    def plot_energy_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the energy density.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_scalar_field`

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

    def plot_load_vector(self, **kwargs):
        """Plot the load vector.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Load vector'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.load_vector, **default_kwargs)

    def plot_current_density(self, **kwargs):
        """Plot the current density.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_scalar_field`

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Current density'}
        default_kwargs.update(kwargs)

        return self.plotter.plot_scalar_field(self.current_density / self.length, **default_kwargs)

    def plot_equilines(self, **kwargs):
        """Plot the equilines of the vector potential. This is in direction of the b field.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`~pyrit.problem.solutions.StaticFieldPlotter.plot_scalar_field`
        """
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis',
                          'title': 'Equilines vector potential', 'colors': 'b'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_equilines(self.vector_potential, **default_kwargs)


class MagneticSolutionCartHarmonic(HarmonicSolution):
    """Solution of a harmonic, two-dimensional magnetic problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemCartHarmonic`.
    """

    def consistency_check(self):
        pass


class MagneticSolutionCartTransient(TransientSolution):
    """Solution of a transient, two-dimensional magnetic problem in Cartesian coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemCartTransient`.
    """

    related_static_solution = MagneticSolutionCartStatic

    def plot_monitor_solution(self, key: str, **kwargs):
        pass

    def consistency_check(self):
        pass
