# coding=utf-8
"""Axisymmetric, magnetic 2D solutions

.. sectionauthor:: Bundschuh
"""
# pylint: disable=used-before-assignment,missing-function-docstring,no-member

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Tuple, Type, Union, Dict, Iterable

import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy.sparse import linalg as splinalg
from scipy import sparse

from pyrit import get_logger
from pyrit import region, material
from pyrit import excitation

from pyrit.solution.Solution import Solution, StaticSolution, HarmonicSolution, TransientSolution, get_field, \
    set_field, get_matrix, set_matrix
from pyrit.solution.FieldPlotter import StaticFieldPlotter, HarmonicFieldPlotter, TransientFieldPlotter

if TYPE_CHECKING:
    from pyrit import shapefunction, mesh

logger = get_logger(__name__)


class VectorConverter:
    r"""Convert vectors for axi meshes of different types.

    For axisymmetric edge functions, there are three different kinds of coefficients for basis functions:

        1. **c-vector**: A vector where the values represent coefficients for the basis functions :math:`\vec{w}_i`.
        2. **r-vector**: A vector where the values represent :math:`\int_\Omega \vec{f}\cdot\vec{w}_i\mathbf{d} V`.
           These appear on the right-hand-side
        3. **f-vector**:  A vector where the values represent the value of the field at a node

    """

    @staticmethod
    def c_to_f_vector(c_vector: np.ndarray, axi_mesh: 'mesh.AxiMesh', nodes_not_on_axis: np.ndarray = None) -> \
            np.ndarray:
        """Convert a c-vector to an f-vector.

        Parameters
        ----------
        c_vector : np.ndarray
        axi_mesh : 'mesh.AxiMesh'
        nodes_not_on_axis : np.ndarray, optional
            An array with the nodes that are not on the axis

        Returns
        -------
        f_vector : np.ndarray
        """
        if nodes_not_on_axis is None:
            nodes_not_on_axis = np.setdiff1d(np.arange(axi_mesh.num_node), axi_mesh.nodes_on_axis)
        f_vector = np.zeros_like(c_vector)
        f_vector[nodes_not_on_axis] = 1 / (2 * np.pi * axi_mesh.node[nodes_not_on_axis, 0]) * c_vector[
            nodes_not_on_axis]
        return f_vector

    @staticmethod
    def f_to_c_vector(f_vector: np.ndarray, axi_mesh: 'mesh.AxiMesh') -> np.ndarray:
        """Convert an f-vector to a c-vector.

        Parameters
        ----------
        f_vector : np.ndarray
        axi_mesh : 'mesh.AxiMesh'

        Returns
        -------
        c_vector : np.ndarray
        """
        return (2 * np.pi * axi_mesh.node[:, 0]) * f_vector

    @staticmethod
    def r_to_c_vector(r_vector: np.ndarray, axi_mesh: 'mesh.AxiMesh', mass_one: sparse.csr_matrix,
                      nodes_not_on_axis: np.ndarray = None) -> np.ndarray:
        """Convert an r-vector to a c-vector.

        Parameters
        ----------
        r_vector : np.ndarray
        axi_mesh : 'mesh.AxiMesh'
        mass_one : sparse.csr_matrix
            A mass matrix with mass 1.
        nodes_not_on_axis : np.ndarray, optional
            An array with the nodes that are not on the axis

        Returns
        -------
        c_vector : np.ndarray
        """
        if nodes_not_on_axis is None:
            nodes_not_on_axis = np.setdiff1d(np.arange(axi_mesh.num_node), axi_mesh.nodes_on_axis)
        if not sparse.isspmatrix_csr(mass_one):
            mass_one = mass_one.tocsr()
        matrix_shrink = mass_one[np.ix_(nodes_not_on_axis, nodes_not_on_axis)]
        s_vector_shrink = r_vector[nodes_not_on_axis]

        tmp = splinalg.spsolve(matrix_shrink, s_vector_shrink)
        c_vector = np.zeros(axi_mesh.num_node, dtype=r_vector.dtype)
        c_vector[nodes_not_on_axis] = tmp
        return c_vector

    @staticmethod
    def c_to_r_vector(c_vector: np.ndarray, mass_one: sparse.csr_matrix) -> np.ndarray:
        """Convert a c-vector to an r-vector.

        Parameters
        ----------
        c_vector : np.ndarray
        mass_one : sparse.csr-matrix

        Returns
        -------
        r_vector : np.ndarray
        """
        return mass_one @ c_vector

    @staticmethod
    def f_to_r_vector(f_vector: np.ndarray, mesh: 'mesh.AxiMesh', mass_one: sparse.csr_matrix) -> np.ndarray:
        """Convert an f-vector to an r-vector.

        Parameters
        ----------
        f_vector : np.ndarray
        mesh : mesh.AxiMesh
        mass_one : sparse.csr_matrix
            A mass matrix with mass 1

        Returns
        -------
        r_vector : np.ndarray
        """
        return VectorConverter.c_to_r_vector(VectorConverter.f_to_c_vector(f_vector, mesh), mass_one)

    @staticmethod
    def r_to_f_vector(r_vector: np.ndarray, axi_mesh: 'mesh.AxiMesh', mass_one: sparse.csr_matrix,
                      nodes_not_on_axis: np.ndarray = None) -> np.ndarray:
        """Convert an r-vector to an f-vector.

        Parameters
        ----------
        r_vector :
        axi_mesh : mesh.AxiMesh
        mass_one : sparse.csr_matrix
            A mass matrix with mass 1.
        nodes_not_on_axis : np.ndarray, optional
            An array with the nodes that are not on the axis

        Returns
        -------
        f_vector : np.ndarray
        """
        return VectorConverter.c_to_f_vector(VectorConverter.r_to_c_vector(r_vector, axi_mesh, mass_one,
                                                                           nodes_not_on_axis),
                                             axi_mesh, nodes_not_on_axis)


# noinspection PyAttributeOutsideInit
class MagneticSolutionAxiBase(Solution, ABC):
    """Abstract class for the solution of static and harmonic, magnetic, axisymmetric solutions."""

    solution_identifier: str = 'Axisymmetric magnetic solution base'

    #: Protected attribute names that are used for properties
    _property_attributes = ['_nodes_not_on_axis', '_vector_potential_field_values', '_b_field', '_h_field',
                            '_curlcurl_matrix', '_mass', '_mass_one', '_load_vector']
    ignore_for_saving = _property_attributes.copy()

    @abstractmethod
    def __init__(self, description: str, vector_potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_edge_shape_function: 'shapefunction.TriAxisymmetricEdgeShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations):
        """Base class for the solution of static and harmonic, magnetic, axisymmetric solutions.

        Parameters
        ----------
        description : str
            A description of the problem.
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`.
        tri_axi_edge_shape_function : TriAxisymmetricEdgeShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricEdgeShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        """
        super().__init__(description, vector_potential, None, None, regions, materials, excitations)
        self.mesh: mesh.AxiMesh = axi_mesh
        self.shape_function: shapefunction.TriAxisymmetricEdgeShapeFunction = tri_axi_edge_shape_function
        self.plotter = None

    def field_vector(self, c_vector: np.ndarray) -> np.ndarray:
        """Wrapper around :py:meth:`VectorConverter.c_to_f_vector` with class specific defaults.

        Parameters
        ----------
        c_vector : np.ndarray

        Returns
        -------
        f_vector : np.ndarray
        """
        return VectorConverter.c_to_f_vector(c_vector, self.mesh, self.nodes_not_on_axis)

    def project_to_c_vector(self, r_vector: np.ndarray) -> np.ndarray:
        """Wrapper around :py:meth:`VectorConverter.r_to_c_vector` with class specific defaults.

        Parameters
        ----------
        r_vector : np.ndarray

        Returns
        -------
        c_vector : np.ndarray
        """
        return VectorConverter.r_to_c_vector(r_vector, self.mesh, self.mass_one, self.nodes_not_on_axis)

    @property
    def nodes_not_on_axis(self) -> np.ndarray:
        """An array with the node indices that are not on the z-axis."""
        if self._nodes_not_on_axis is None:
            self._nodes_not_on_axis = np.setdiff1d(np.arange(self.mesh.num_node), self.mesh.nodes_on_axis)
        return self._nodes_not_on_axis

    @nodes_not_on_axis.setter
    def nodes_not_on_axis(self, nodes_not_on_axis: np.ndarray):
        self._nodes_not_on_axis = nodes_not_on_axis

    @property
    def vector_potential(self) -> np.ndarray:
        """The vector potential. These are the coefficients for the basis functions."""
        return self._solution

    @vector_potential.setter
    def vector_potential(self, vector_potential: np.ndarray):
        self._solution = vector_potential

    @property
    def vector_potential_field_values(self):
        """The field values of the vector potential."""
        if self._vector_potential_field_values is None:
            self._vector_potential_field_values = self.field_vector(self.vector_potential)
        return self._vector_potential_field_values

    @vector_potential_field_values.setter
    def vector_potential_field_values(self, vector_potential_field_values):
        self._vector_potential_field_values = vector_potential_field_values

    @property
    def b_field(self) -> np.ndarray:
        """Magnetic flux density. Array with two values (r- and z-component) per mesh element."""
        if self._b_field is None:
            self._b_field = self.shape_function.curl(self.vector_potential)
        return self._b_field

    @b_field.setter
    def b_field(self, b_field: np.ndarray):
        self._b_field = b_field

    @property
    def h_field(self) -> np.ndarray:
        """Magnetic field strength. Array with two values (r- and z-component) per mesh element."""
        if self._h_field is None:
            reluctivity_per_element = self.mesh.value_per_element(self.regions, self.materials, material.Reluctivity)
            self._h_field = np.tile(np.atleast_2d(reluctivity_per_element).T, 2) * self.b_field
        return self._h_field

    @h_field.setter
    def h_field(self, h_field: np.ndarray):
        self._h_field = h_field

    @property
    def curlcurl_matrix(self) -> sparse.coo_matrix:
        """The curl-curl matrix."""
        if self._curlcurl_matrix is None:
            self._curlcurl_matrix = self.shape_function.curlcurl_operator(self.regions, self.materials,
                                                                          material.Reluctivity)
        return self._curlcurl_matrix

    @curlcurl_matrix.setter
    def curlcurl_matrix(self, curlcurl_matrix: sparse.coo_matrix):
        self._curlcurl_matrix = curlcurl_matrix

    @property
    def mass(self) -> sparse.coo_matrix:
        """The mass matrix from the conductivity"""
        if self._mass is None:
            self._mass = self.shape_function.mass_matrix(self.regions, self.materials, material.Conductivity)
        return self._mass

    @mass.setter
    def mass(self, mass: sparse.coo_matrix):
        self._mass = mass

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
    def mass_one(self) -> sparse.csr_matrix:
        """The mass matrix with mass 1."""
        if self._mass_one is None:
            self._mass_one = self.shape_function.mass_matrix(1).tocsr()
        return self._mass_one

    @mass_one.setter
    def mass_one(self, mass_one: sparse.csr_matrix):
        self._mass_one = mass_one

    @property
    def current_density(self) -> np.ndarray:
        """The current density."""
        return self.project_to_c_vector(self.load_vector)

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
            The plot_type of the plot. See :py:func:`plot_field`
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
        if plot_type == 'line':
            kwargs.setdefault('colors', 'blue')
            kwargs.setdefault('num_lines', 20)
            return self.plotter.plot_equilines(self.vector_potential / (2 * np.pi), **kwargs)
        return self.plotter.plot_vector_field(self.b_field, plot_type, **kwargs)

    def plot_h_field(self, plot_type: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the magnetic field strength.

        Parameters
        ----------
        plot_type : str
            The plot_type of the plot. See :py:func:`plot_field`
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

    def plot_load_vector(self, **kwargs):
        """Plot the load vector.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Load vector'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.load_vector, **default_kwargs)

    def plot_current_density(self, **kwargs):
        """Plot the current density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Current density'}
        default_kwargs.update(kwargs)

        return self.plotter.plot_scalar_field(self.field_vector(self.current_density), **default_kwargs)

    def plot_equilines(self, **kwargs):
        """Plot the equilines of the vector potential. This is in direction of the b field."""
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis',
                          'title': 'Equilines vector potential', 'colors': 'b'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_equilines(self.vector_potential / (2 * np.pi), **default_kwargs)


# noinspection PyAttributeOutsideInit
class MagneticSolutionAxiStatic(StaticSolution, MagneticSolutionAxiBase):
    """Solution of a static, two-dimensional magnetic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemAxiStatic`.
    """

    solution_identifier: str = 'Static, axisymmetric magnetic solution'

    def __init__(self, description: str, vector_potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_edge_shape_function: 'shapefunction.TriAxisymmetricEdgeShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations):
        super().__init__(description, vector_potential, axi_mesh, tri_axi_edge_shape_function, regions, materials,
                         excitations)
        self.plotter = StaticFieldPlotter(axi_mesh)

        self.consistency_check()

    def consistency_check(self):
        if not self.vector_potential.shape == (self.mesh.num_node,):
            logger.critical("The vector potential has a unexpected format. It can't be used to compute other data.")
        if not issubclass(self.vector_potential.dtype.type, (np.integer, np.floating)):
            logger.warning("The dtype of the vector potential is '%s'. This is not expected.",
                           self.vector_potential.dtype)

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

        mean_radii = np.mean(self.mesh.node[self.mesh.elem2node[indices], 0], axis=1)
        energy = 2 * np.pi * np.sum(self.energy_density[indices] * self.mesh.elem_area[indices] * mean_radii)

        return energy

    def plot_energy_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the energy density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Energy density'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.energy_density, **default_kwargs)


# noinspection PyAttributeOutsideInit
class MagneticSolutionAxiHarmonic(HarmonicSolution, MagneticSolutionAxiBase):
    """Solution of a harmonic, two-dimensional magnetic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemAxiHarmonic`.
    """

    solution_identifier: str = 'Harmonic, axisymmetric magnetic solution'

    def __init__(self, description: str, vector_potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_edge_shape_function: 'shapefunction.TriAxisymmetricEdgeShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations, frequency: float):
        """A harmonic, magnetic, axisymmetric solution.

        Parameters
        ----------
        description : str
            A description of the problem.
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`.
        tri_axi_edge_shape_function : TriAxisymmetricEdgeShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricEdgeShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        frequency : float
            The frequency of the solution.
        """
        super().__init__(description, vector_potential, axi_mesh, tri_axi_edge_shape_function, regions, materials,
                         excitations, frequency)
        self.plotter = HarmonicFieldPlotter(axi_mesh, frequency)

        self.consistency_check()

    def consistency_check(self):
        if not self.vector_potential.shape == (self.mesh.num_node,):
            logger.critical("The vector potential has a unexpected format. It can't be used to compute other data.")
        if not issubclass(self.vector_potential.dtype.type, np.complexfloating):
            logger.warning("The dtype of the vector potential is '%s'. This is expected to be complex.",
                           self.vector_potential.dtype)
        if self.frequency < 0:
            logger.error("The frequency is negative.")
        if self.frequency == 0:
            logger.warning("The frequency is zero.")

    def energy(self, **kwargs) -> float:
        """The energy in the whole domain.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`get_field_at_phase`. Specify phase or time.

        Returns
        -------
        energy : float
            The energy in the whole domain.
        """
        vec_pot = self.get_field_at_phase(self.frequency, self.vector_potential, **kwargs)[0]
        return 0.5 * vec_pot @ self.curlcurl_matrix @ vec_pot

    def energy_density(self, **kwargs) -> np.ndarray:
        """The energy density. Array with one value per mesh node.

        Parameters
        ----------
        kwargs :
            Passed to :py:meth:`~HarmonicFieldPlotter.get_field_at_phase`. Specify phase or time.

        Returns
        -------
        energy_density : np.ndarray
            An array with the energy density per element.
        """
        b_field = self.get_field_at_phase(self.frequency, self.b_field, **kwargs)[0]
        h_field = self.get_field_at_phase(self.frequency, self.h_field, **kwargs)[0]
        return 0.5 * np.sum(b_field * h_field, axis=1)

    def energy_in_regions(self, *region_ids: int, **kwargs) -> float:
        """Compute the energy in certain regions.

        Parameters
        ----------
        region_ids : int
            Region IDs of the regions to compute the energy in.

        Returns
        -------
        energy : float
            The energy
        kwargs :
            Passed to :py:meth:`energy_density`. Specify phase or time.

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

        mean_radii = np.mean(self.mesh.node[self.mesh.elem2node[indices], 0], axis=1)
        energy = 2 * np.pi * np.sum(self.energy_density(**kwargs)[indices] * self.mesh.elem_area[indices] * mean_radii)

        return energy

    @property
    def source_current_density(self) -> np.ndarray:
        """The source current density."""
        return VectorConverter.r_to_c_vector(self.load_vector, self.mesh, self.mass_one, self.nodes_not_on_axis)

    @property
    def eddy_current_density(self):
        """The eddy-current density."""
        return - 1j * self.omega * VectorConverter.r_to_c_vector(self.mass @ self.vector_potential, self.mesh,
                                                                 self.mass_one, self.nodes_not_on_axis)

    @property
    def current_density(self) -> np.ndarray:
        """The whole current density."""
        return self.source_current_density + self.eddy_current_density

    def plot_source_current_density(self, **kwargs):
        """Plot the source current density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Source current density'}
        default_kwargs.update(kwargs)

        return self.plotter.plot_scalar_field(self.field_vector(self.source_current_density), **default_kwargs)

    def plot_eddy_current_density(self, **kwargs):
        """Plot the eddy current density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Eddy current density'}
        default_kwargs.update(kwargs)

        return self.plotter.plot_scalar_field(self.field_vector(self.eddy_current_density), **default_kwargs)

    def plot_energy_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the energy density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Energy density'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.energy_density(**kwargs), **default_kwargs)


# noinspection PyAttributeOutsideInit,PyUnresolvedReferences,PyUnusedLocal
class MagneticSolutionAxiTransient(TransientSolution):
    """Solution of a transient, two-dimensional magnetic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.MagneticProblemAxiTransient`.

    For all get- and set-methods with time and index in the parameters applies the following: If both are not given (or
    None) get or set the whole value. If only one in not None, get or set the value at the specific time instance. If
    both are not None, an error is raised because of ambiguity.
    """

    solution_identifier: str = 'Transient, axisymmetric magnetic solution'
    related_static_solution = MagneticSolutionAxiStatic

    _property_attributes = ['_nodes_not_on_axis', '_mass_one']

    _derived_attributes = ['_vector_potential_field_values', '_b_field', '_h_field', '_energy', '_energy_density',
                           '_energy_in_regions', '_eddy_current_density', '_source_current_density', '_current_density']

    _matrix_attributes = ['_curlcurl_matrix', '_mass_matrix', '_mass_matrix_one', '_load_vector']

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, vector_potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_edge_shape_function: 'shapefunction.TriAxisymmetricEdgeShapeFunction', regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations, time_steps: np.ndarray,
                 monitor_solutions: Dict[str, Tuple[np.ndarray, np.ndarray]]):
        """Constructor.

        Parameters
        ----------
        description :
        vector_potential :
        axi_mesh :
        tri_axi_edge_shape_function :
        regions :
        materials :
        excitations :
        """
        super().__init__(description, vector_potential, axi_mesh, tri_axi_edge_shape_function, regions, materials,
                         excitations, time_steps, monitor_solutions)
        self.mesh: 'mesh.AxiMesh' = axi_mesh
        self.shape_function: 'shapefunction.TriAxisymmetricEdgeShapeFunction' = tri_axi_edge_shape_function
        self.plotter = TransientFieldPlotter(axi_mesh)

        self.consistency_check()

    def consistency_check(self):
        if not self.vector_potential().shape == (len(self.time_steps), self.mesh.num_node):
            if self.vector_potential().shape == (self.mesh.num_node, len(self.time_steps)):
                logger.warning("You passed the vector potential in wrong format. It is going to be transformed.")
                self.set_vector_potential(self.vector_potential().transpose())
            else:
                logger.critical("The vector potential has a unexpected format. It can't be used to compute other data.")
        if not issubclass(self.vector_potential().dtype.type, (np.integer, np.floating)):
            logger.warning("The dtype of the vector potential is '%s'. This is not expected.",
                           self.vector_potential().dtype)
        if len(self.time_steps.shape) != 1:
            logger.critical("The shape of 'time_steps' is not as expected.")

    @property
    def nodes_not_on_axis(self):
        """An array with the node indices that are not on the z-axis."""
        if self._nodes_not_on_axis is None:
            self._nodes_not_on_axis = np.setdiff1d(np.arange(self.mesh.num_node), self.mesh.nodes_on_axis)
        return self._nodes_not_on_axis

    @nodes_not_on_axis.setter
    def nodes_not_on_axis(self, nodes_not_on_axis):
        self._nodes_not_on_axis = nodes_not_on_axis

    def field_vector(self, c_vector):
        """Wrapper around :py:meth:`VectorConverter.c_to_f_vector` with class specific defaults.

        Parameters
        ----------
        c_vector : np.ndarray

        Returns
        -------
        f_vector : np.ndarray
        """
        return VectorConverter.c_to_f_vector(c_vector, self.mesh, self.nodes_not_on_axis)

    def project_to_c_vector(self, r_vector):
        """Wrapper around :py:meth:`VectorConverter.r_to_c_vector` with class specific defaults.

        Parameters
        ----------
        r_vector : np.ndarray

        Returns
        -------
        c_vector : np.ndarray
        """
        return VectorConverter.r_to_c_vector(r_vector, self.mesh, self.mass_one, self.nodes_not_on_axis)

    def vector_potential(self, *, time=None, index=None) -> np.ndarray:
        """The vector potential. These are the coefficients for the basis functions.

        Parameters
        ----------
        time : float
            A time instance.
        index : int
            The index of a time instance.

        Returns
        -------
        vector_potential : np.ndarray
        """
        return self._get_solution(time=time, index=index)

    def set_vector_potential(self, vector_potential: np.ndarray, *, time=None, index=None):
        """Set the vector potential.

        Parameters
        ----------
        vector_potential : np.ndarray
            The vector potential to set
        time : float
            A time instance.
        index : int
            The index of a time instance.
        """
        self._set_solution(vector_potential, time=time, index=index)

    # region Derived values

    @get_field('_vector_potential_field_values', 'vector potential', 'ndarray')
    def vector_potential_field_values(self, *, time=None, index=None):
        return self.field_vector(self.vector_potential(index=index))

    @set_field('_vector_potential_field_values', 'vector potential', 'ndarray')
    def set_vector_potential_field_values(self, vector_potential_field_values: Union[np.ndarray, Dict[int, np.ndarray]],
                                          *, time=None, index=None):
        self._check_set_types(vector_potential_field_values, np.ndarray)

    @get_field('_b_field', 'magnetic flux density', 'ndarray')
    def b_field(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        return self.shape_function.curl(self.vector_potential(index=index))

    @set_field('_b_field', 'magnetic flux density', 'ndarray')
    def set_b_field(self, b_field: Union[np.ndarray, Dict[int, np.ndarray]], *, time=None, index=None):
        self._check_set_types(b_field, np.ndarray)

    @get_field('_h_field', 'magnetic field strength', 'ndarray')
    def h_field(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        reluctivity_per_element = self.mesh.value_per_element(self.regions, self.materials, material.Reluctivity)
        return np.tile(np.atleast_2d(reluctivity_per_element).T, 2) * self.b_field(index=index)

    @set_field('_h_field', 'magnetic field strength', 'ndarray')
    def set_h_field(self, h_field: Union[np.ndarray, Dict[int, np.ndarray]], *, time=None, index=None):
        self._check_set_types(h_field, np.ndarray)

    @get_field('_energy', 'magnetic energy', 'float')
    def energy(self, *, time=None, index=None) -> Union[float, Dict[int, float]]:
        vector_potential = self.vector_potential(index=index)
        curlcurl_matrix = self.curlcurl_matrix(index=index)
        return 0.5 * vector_potential @ curlcurl_matrix @ vector_potential

    @set_field('_energy', 'magnetic energy', 'float')
    def set_energy(self, energy: Union[float, Dict[int, float]], *, time=None, index=None):
        self._check_set_types(energy, float)

    @get_field('_energy_density', 'magnetic energy density', 'ndarray')
    def energy_density(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        b_field = self.b_field(index=index)
        h_field = self.h_field(index=index)
        return 0.5 * np.sum(b_field * h_field, axis=1)

    @set_field('_energy_density', 'magnetic energy density', 'ndarray')
    def set_energy_density(self, energy_density: Union[np.ndarray, Dict[int, np.ndarray]], *, time=None, index=None):
        self._check_set_types(energy_density, np.ndarray)

    @get_field('_eddy_current_density', 'eddy current density', 'ndarray')
    def eddy_current_density(self, *, time=None, index=None):
        if index == 0:
            return np.zeros(self.mesh.num_node)

        delta_t = self.time_steps[index] - self.time_steps[index - 1]
        delt_a = (self.vector_potential(index=index) - self.vector_potential(index=index - 1)) / delta_t
        mass = self.mass_matrix(index=index)
        return -1 * self.project_to_c_vector(mass @ delt_a)

    @set_field('_eddy_current_density', 'eddy current density', 'ndarray')
    def set_eddy_current_density(self, eddy_current_density: Union[np.ndarray, Dict[int, np.ndarray]], *, time=None,
                                 index=None):
        self._check_set_types(eddy_current_density, np.ndarray)

    @get_field('_source_current_density', 'source current density', 'ndarray')
    def source_current_density(self, *, time=None, index=None) -> np.ndarray:
        return self.project_to_c_vector(self.load_vector(time=time, index=index).toarray())

    @get_field('_current_density', 'current density', 'ndarray')
    def current_density(self, *, time=None, index=None):
        index = self.get_index(time, index)
        return self.source_current_density(index=index) + self.eddy_current_density(index=index)

    def energy_in_regions(self, *region_ids: int, time=None, index=None) -> float:
        """Compute the energy in certain regions.

        Parameters
        ----------
        region_ids : int
            Region IDs of the regions to compute the energy in.
        time : float, optional
            Time instance.
        index : int, optional
            The index of a time instance.

        Returns
        -------
        energy : float
            The energy

        Raises
        ------
        ValueError : When the region of one ID is not two-dimensional.
        """
        if time is None and index is None:
            return self._energy_in_regions
        index = self.get_index(time, index)

        # If the time index has not calculated yet, add an empty dict
        if index not in self._energy_in_regions.keys():
            self._energy_in_regions[index] = {}

        indices_per_region = []  # List of indices of the elements in the regions
        for region_id in region_ids:
            if not self.regions.get_regi(region_id).dim == 2:
                raise ValueError(f"The region '{region_id}' does not belong to a region of the expected dimension.")
            if region_id not in self._energy_in_regions[index].keys():
                indices = indices_per_region.append(np.where(self.mesh.elem2regi == region_id)[0])
                self._energy_in_regions[index][region_id] = 2 * np.pi * np.sum(
                    self.energy_density(index=index)[indices] * self.mesh.elem_area[indices] * self.mesh.node[0, :])

        return sum((self._energy_in_regions[index][region_id] for region_id in region_ids))

    # endregion

    # region Matrix and vector

    @get_matrix('_curlcurl_matrix', 'curl-curl matrix')
    def curlcurl_matrix(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.Reluctivity)
        time_dependent = self.materials.is_time_dependent(material.Reluctivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.Reluctivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.Reluctivity)

        return lin and not time_dependent, lambda: self.shape_function.curlcurl_operator(self.regions, self.materials,
                                                                                         material.Reluctivity)

    @set_matrix('_curlcurl_matrix', 'curl-curl matrix')
    def set_curlcurl_matrix(self, curlcurl_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_mass_matrix', 'mass matrix')
    def mass_matrix(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.Conductivity)
        time_dependent = self.materials.is_time_dependent(material.Conductivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.Conductivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.Conductivity)

        return lin and not time_dependent, lambda: self.shape_function.mass_matrix(self.regions, self.materials,
                                                                                   material.Conductivity)

    @set_matrix('_mass_matrix', 'mass matrix')
    def set_mass_matrix(self, mass_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_load_vector', 'load vector')
    def load_vector(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.excitations.is_linear
        time_dependent = self.excitations.is_time_dependent

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.excitations.update('solution', static_solution)
        if time_dependent:
            self.excitations.update('time', self.time_steps[index])

        return lin and not time_dependent, lambda: self.shape_function.load_vector(self.regions, self.excitations)

    @set_matrix('_load_vector', 'load vector')
    def set_load_vector(self, load_vector: sparse.coo_matrix, *, time=None, index=None):
        pass

    @property
    def mass_one(self) -> sparse.csr_matrix:
        """The mass matrix with mass 1."""
        if self._mass_one is None:
            self._mass_one = self.shape_function.mass_matrix(1).tocsr()
        return self._mass_one

    @mass_one.setter
    def mass_one(self, mass_one: sparse.csr_matrix):
        self._mass_one = mass_one

    # endregion

    # region Plot on time instance and animate

    def plot_vector_potential(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the vector potential.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Vector potential'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.vector_potential_field_values(index=index),
                                              self.time_steps[index], **kwargs)

    def animate_vector_potential(self, *, times: Iterable[float] = None,
                                 indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the vector potential.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Vector potential'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.animate_scalar_field((self.vector_potential_field_values(index=i) for i in indices),
                                                 **kwargs)

    def plot_b_field(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the magnetic flux density.

        Parameters
        ----------
        plot_type : str
            The plot_type of the plot. See :py:func:`plot_field`
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to the plot method.

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_vector_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic flux density $\vec{B}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_vector_field(self.b_field(index=index), plot_type, self.time_steps[index], **kwargs)

    def animate_b_field(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                        indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the magnetic flux density.

        Parameters
        ----------
        plot_type : str
            The plot_type of the vector plot. See :py:meth:`TransientFieldPlotter.animate_vector_field`.
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_vector_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_vector_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic flux density $\vec{B}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.b_field(index=i) for i in indices), plot_type, **kwargs)

    def plot_h_field(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the magnetic field strength.

        Parameters
        ----------
        plot_type : str
            The plot_type of the plot. See :py:func:`plot_field`
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to the plot method.

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_vector_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic field strength $\vec{H}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.h_field(index=index), plot_type, self.time_steps[index], **kwargs)

    def animate_h_field(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                        indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the magnetic field strength.

        Parameters
        ----------
        plot_type : str
            The plot_type of the vector plot. See :py:meth:`TransientFieldPlotter.animate_vector_field`.
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_vector_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_vector_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic flux density $\vec{B}$"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.h_field(index=i) for i in indices), plot_type, **kwargs)

    def plot_energy_density(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the energy density.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Energy density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.energy_density(index=index), self.time_steps[index], **kwargs)

    def animate_energy_density(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                               **kwargs):
        """Animate the energy_density.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Magnetic energy density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.energy_density(index=i) for i in indices), **kwargs)

    def plot_load_vector(self, *, time=None, index=None, **kwargs):
        """Plot the load vector.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Load vector'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.load_vector(index=index).toarray(), self.time_steps[index], **kwargs)

    def animate_load_vector(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                            **kwargs):
        """Animate the load vector.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        if self.excitations.is_linear and not self.excitations.is_time_dependent:
            logger.warning("The load vector is linear and independent of time. An animation does not make sense.")
        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': r"Load vector"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.load_vector(index=i).toarray() for i in indices), **kwargs)

    def plot_source_current_density(self, *, time=None, index=None, **kwargs):
        """Plot the source current density.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Source current density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_scalar_field(self.field_vector(self.source_current_density(index=index)),
                                              self.time_steps[index], **kwargs)

    def animate_source_current_density(self, *, times: Iterable[float] = None,
                                       indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the source current density.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': "Source current density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.field_vector(self.source_current_density(index=i)) for i in indices),
                                          **kwargs)

    def plot_eddy_current_density(self, *, time=None, index=None, **kwargs):
        """Plot the eddy current density.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Eddy current density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_scalar_field(self.field_vector(self.eddy_current_density(index=index)),
                                              self.time_steps[index], **kwargs)

    def animate_eddy_current_density(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                                     **kwargs):
        """Animate the eddy source current density.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': "Eddy current density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.field_vector(self.eddy_current_density(index=i)) for i in indices),
                                          **kwargs)

    def plot_current_density(self, *, time=None, index=None, **kwargs):
        """Plot the current density.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_scalar_field`
        """
        index = self.get_index(time, index)

        default_kwargs = {'plot_3d': False, 'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Current density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_scalar_field(self.field_vector(self.current_density(index=index)),
                                              self.time_steps[index], **kwargs)

    def animate_current_density(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                                **kwargs):
        """Animate the current density.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': "Current density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.field_vector(self.current_density(index=i)) for i in indices), **kwargs)

    def plot_equilines(self, *, time=None, index=None, **kwargs):
        """Plot the equilines of the vector potential. This is in direction of the b field.

        Parameters
        ----------
        time : float, optional
            A time instance.
        index : int, optional
            The index of a time instance.
        kwargs :
            Passed to :py:meth:`TriMesh.plot_scalar_field`

        Returns
        -------
        Return of :py:meth:`~TransientFieldPlotter.plot_equilines`
        """
        index = self.get_index(time, index)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': 'Equilines vector potential',
                          'colors': 'b'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_equilines(self.vector_potential_field_values(index=index),
                                           self.time_steps[index], **kwargs)

    def animate_equilines(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the equilines of the vector potential.

        Parameters
        ----------
        times : Iterable[float]
            An iterable with time instances. See :py:func:`get_animation_indices`.
        indices : Union[Iterable[int], int]
            An iterable with indices. See :py:func:`get_animation_indices`.
        kwargs :
            Keyword arguments passed to :py:meth:`TransientFieldPlotter.animate_scalar_field`.

        Returns
        -------
        Return of :py:meth:`TransientFieldPlotter.animate_scalar_field`
        """
        indices = self.get_animation_indices(times, indices)

        default_kwargs = {'x_label': 'x axis', 'y_label': 'y axis', 'title': "Equilines vector potential"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.vector_potential_field_values(index=i) for i in indices), **kwargs)

    # endregion

    def plot_monitor_solution(self, key: str, **kwargs):
        if key not in self.monitor_solutions:
            raise KeyError(f"Key '{key}' not known in monitor solutions")
        return self.plotter.plot_over_time(self.monitor_solutions[key][0], self.monitor_solutions[key][1], **kwargs)
