# coding=utf-8
"""Current flow solution in axisymmetric coordinates

.. sectionauthor:: Bundschuh, Ruppert
"""
# pylint: disable=used-before-assignment,missing-function-docstring,no-member,attribute-defined-outside-init

from typing import TYPE_CHECKING, Tuple, Type, Union, Dict, Iterable

import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy import sparse

from pyrit import get_logger
from pyrit import region, material
from pyrit import excitation

from pyrit.solution.FieldPlotter import StaticFieldPlotter, TransientFieldPlotter
from pyrit.solution.Solution import StaticSolution, HarmonicSolution, TransientSolution, get_field, set_field, \
    get_matrix, set_matrix

if TYPE_CHECKING:
    from pyrit import shapefunction, mesh

logger = get_logger(__name__)

__all__ = ['CurrentFlowSolutionAxiStatic', 'CurrentFlowSolutionAxiHarmonic', 'CurrentFlowSolutionAxiTransient']


class CurrentFlowSolutionAxiStatic(StaticSolution):
    """The solution of a stationary current problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemAxiStatic`.
    """

    #: Protected attribute names that are used for properties
    _derived_attributes = ['_e_field', '_e_field_abs',
                           '_joule_loss_power_density', '_current_density',
                           '_current_density_abs', '_conductivity_per_element',
                           '_joule_loss_power', '_joule_loss_power_in_regions', '_conductivity_per_element']
    _matrix_attributes = ['_divgrad_matrix_sigma', '_divgrad_matrix_sigmad',
                          '_load_vector_electric']

    _property_attributes = _derived_attributes.copy() + _matrix_attributes.copy()

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_nodal_shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction',
                 regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations):
        """Base class for the solution of static and harmonic, magnetic, axisymmetric solutions.

        Parameters
        ----------
        description : str
            A description of the problem.
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`.
        tri_axi_nodal_shape_function : TriAxisymmetricNodalShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricEdgeShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        """
        super().__init__(description, potential, None, None, regions, materials, excitations)
        self.mesh: mesh.AxiMesh = axi_mesh
        self.shape_function: shapefunction.TriAxisymmetricNodalShapeFunction = tri_axi_nodal_shape_function
        self.plotter = StaticFieldPlotter(axi_mesh)
        self.consistency_check()
        self.update_material_electric()

    def update_material_electric(self):
        # Write solution to nonlinear materials and excitations to ensure consistency.
        self.materials.update(name="solution", val=self, prop_class=material.Conductivity)
        self.materials.update(name="solution", val=self, prop_class=material.DifferentialConductivity)
        self.materials.update(name="potential", val=self.potential, prop_class=material.ThermalConductivity)
        self.materials.update(name="solution", val=self, prop_class=material.Permittivity)
        self.materials.update(name="solution", val=self, prop_class=material.DifferentialPermittivity)
        self.materials.update(name="potential", val=self.potential, prop_class=material.VolumetricHeatCapacity)
        self.excitations.update(name="solution", val=self)

    def ensure_consistency_potential(self):
        # Function that ensures consistency by resetting all values derived from the potential.
        for prop_name in CurrentFlowSolutionAxiStatic._derived_attributes:
            setattr(self, prop_name, None)
        self.update_material_electric()

        if not self.materials.is_linear(material.Conductivity):
            self._divgrad_matrix_sigma = None
        if not self.materials.is_linear(material.DifferentialConductivity):
            self._divgrad_matrix_sigmad = None
        if not self.excitations.is_linear:
            self._load_vector_electric = None

    @property
    def potential(self) -> np.ndarray:
        """The electric potential. These are the coefficients for the basis functions."""
        return self._solution

    @potential.setter
    def potential(self, potential: np.ndarray):
        self._solution = potential
        self.ensure_consistency_potential()

    @property
    def e_field(self) -> np.ndarray:
        """The electric field. Array with two values per mesh element."""
        if self._e_field is None:
            self._e_field = - self.shape_function.gradient(self.potential)
        return self._e_field

    @e_field.setter
    def e_field(self, e_field: np.ndarray):
        self._e_field = e_field

    @property
    def e_field_abs(self) -> np.ndarray:
        """Absolute value of the electric field. Array with one value per mesh element."""
        if self._e_field_abs is None:
            self._e_field_abs = np.linalg.norm(self.e_field, axis=1)
        return self._e_field_abs

    @e_field_abs.setter
    def e_field_abs(self, e_field_abs: np.ndarray):
        self._e_field_abs = e_field_abs

    @property
    def conductivity_per_element(self) -> np.ndarray:
        """Electric conductivity per element. Array with one value per mesh element."""
        if self._conductivity_per_element is None:
            self._conductivity_per_element = self.mesh.value_per_element(self.regions, self.materials,
                                                                         prop_class=material.Conductivity)

        return self._conductivity_per_element

    @conductivity_per_element.setter
    def conductivity_per_element(self, sigma: np.ndarray):
        self._conductivity_per_element = sigma

    @property
    def current_density(self) -> np.ndarray:
        """Electric current density. Array with two values per mesh element."""
        if self._current_density is None:
            self._current_density = self.conductivity_per_element[:, None] * self.e_field
        return self._current_density

    @current_density.setter
    def current_density(self, current_density):
        self._current_density = current_density

    @property
    def current_density_abs(self) -> np.ndarray:
        if self._current_density_abs is None:
            self._current_density_abs = self.conductivity_per_element * self.e_field_abs
        return self._current_density_abs

    @current_density_abs.setter
    def current_density_abs(self, current_density_abs):
        """Absolute value of the current density. Array with one value per mesh element."""
        self._current_density_abs = current_density_abs

    @property
    def joule_loss_power_density(self) -> np.ndarray:
        """Joule loss power density. Array with one value per mesh element."""
        if self._joule_loss_power_density is None:
            self._joule_loss_power_density = np.sum(self.current_density * self.e_field, axis=1)
        return self._joule_loss_power_density

    @joule_loss_power_density.setter
    def joule_loss_power_density(self, pld):
        self._joule_loss_power_density = pld

    @property
    def load_vector_electric(self) -> sparse.coo_matrix:
        """The load vector."""
        if self._load_vector_electric is None:
            self._load_vector_electric = self.shape_function.load_vector(self.regions, self.excitations)
        return self._load_vector_electric

    @load_vector_electric.setter
    def load_vector_electric(self, load: sparse.coo_matrix):
        self._load_vector_electric = load

    @property
    def divgrad_matrix_sigma(self) -> sparse.coo_matrix:
        """The divgrad matrix for the conductivity."""
        if self._divgrad_matrix_sigma is None:
            self._divgrad_matrix_sigma = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                              material.Conductivity)
        return self._divgrad_matrix_sigma

    @divgrad_matrix_sigma.setter
    def divgrad_matrix_sigma(self, divgrad: sparse.coo_matrix):
        self._divgrad_matrix_sigma = divgrad

    @property
    def divgrad_matrix_sigmad(self) -> sparse.coo_matrix:
        """The divgrad matrix for the differential conductivity."""
        if self._divgrad_matrix_sigmad is None:
            self._divgrad_matrix_sigmad = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                               material.DifferentialConductivity)
        return self._divgrad_matrix_sigmad

    @divgrad_matrix_sigmad.setter
    def divgrad_matrix_sigmad(self, divgrad: sparse.coo_matrix):
        self._divgrad_matrix_sigmad = divgrad

    @property
    def joule_loss_power(self):
        """The joule loss power in the whole domain."""
        if self._joule_loss_power is None:
            self._joule_loss_power = self.potential @ self.divgrad_matrix_sigma @ self.potential
        return self._joule_loss_power

    @joule_loss_power.setter
    def joule_loss_power(self, joule_loss_power):
        self._joule_loss_power = joule_loss_power

    def joule_loss_power_in_regions(self, *region_ids: int) -> float:
        """Compute the joule loss power in certain regions.

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
        energy = 2 * np.pi * np.sum(self.joule_loss_power_density[indices] * self.mesh.elem_area[indices] * mean_radii)

        return energy

    def consistency_check(self):
        if self.potential.shape != (self.mesh.num_node,):
            logger.critical("The electric potential has a unexpected format. It can't be used to compute other data.")

    def plot_potential(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the electric potential.

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
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis',
                          'title': 'electric potential in V'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.potential, **kwargs)

    def plot_e_field(self, kind: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the electric field.

        Parameters
        ----------
        kind : str
            The kind of the plot. See :py:func:`plot_field`
        kwargs :
            Passed to the plot method.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"electric field $\vec{E}$ in V/m"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        if kind == 'line':
            kwargs.setdefault('colors', 'blue')
            kwargs.setdefault('num_lines', 20)
            return self.plotter.plot_equilines(self.potential, **kwargs)
        return self.plotter.plot_vector_field(self.e_field, kind, **kwargs)

    def plot_current_density(self, kind: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the electric current density.

        Parameters
        ----------
        kind : str
            The kind of the plot. See :py:func:`plot_field`
        kwargs :
            Passed to the plot method.

        Returns
        -------
        fig: Figure
            A figure object.
        ax: Axes
            A axes object.
        """
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"current density $\vec{J}$ in A/m^2"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.current_density, kind, **kwargs)

    def plot_joule_loss_power_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the joule loss power density.

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
        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis',
                          'title': 'power loss density W/m^3'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.joule_loss_power_density, **default_kwargs)


class CurrentFlowSolutionAxiIntermediate(CurrentFlowSolutionAxiStatic):
    """An intermediate solution of a current flow problem in axisymmetric coordinates.

    An intermediate solution of a current flow problem. It has all properties of a stationary current solution, but
    is extended by properties concerning time-dependencies and the permittivity. It is used for internal
    auxiliary variables in the solve-method of :py:class:`~pyrit.problem.problems.CurrentFlowProblemAxiTransient`.
    """

    #: Protected attribute names that are used for properties
    _derived_attributes = ['_e_field', '_e_field_abs',
                           '_joule_loss_power_density', '_current_density',
                           '_current_density_abs', '_conductivity_per_element',
                           '_joule_loss_power', '_joule_loss_power_in_regions', '_conductivity_per_element']
    _matrix_attributes = ['_divgrad_matrix_sigma', '_divgrad_matrix_sigmad', '_divgrad_matrix_eps',
                          '_divgrad_matrix_epsd', '_load_vector_electric']

    _property_attributes = ['_time'] + _derived_attributes.copy() + _matrix_attributes.copy()

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_nodal_shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction',
                 regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations):
        """Base class for the solution of static and harmonic, magnetic, axisymmetric solutions.

        Parameters
        ----------
        description : str
            A description of the problem.
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`.
        tri_axi_nodal_shape_function : TriAxisymmetricNodalShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricEdgeShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        """
        super().__init__(description, potential, axi_mesh, tri_axi_nodal_shape_function,
                         regions, materials, excitations)

    def ensure_consistency_potential(self):
        # Function that ensures consistency by resetting all values derived from the potential.
        super().ensure_consistency_potential()

        if not self.materials.is_linear(material.Permittivity):
            self._divgrad_matrix_eps = None
        if not self.materials.is_linear(material.DifferentialPermittivity):
            self._divgrad_matrix_epsd = None

    def ensure_consistency_time(self, time):
        # Write time to time-dependent materials and excitations to ensure consistency.
        if self.materials.is_time_dependent(prop_class=material.Permittivity):
            self.materials.update(name="time", val=time, prop_class=material.Permittivity)
            self._divgrad_matrix_eps = None
        if self.materials.is_time_dependent(prop_class=material.DifferentialPermittivity):
            self.materials.update(name="time", val=time, prop_class=material.DifferentialPermittivity)
            self._divgrad_matrix_epsd = None
        if self.materials.is_time_dependent(prop_class=material.Conductivity):
            self.materials.update(name="time", val=time, prop_class=material.Conductivity)
            self._divgrad_matrix_sigma = None
        if self.materials.is_time_dependent(prop_class=material.DifferentialConductivity):
            self.materials.update(name="time", val=time, prop_class=material.DifferentialConductivity)
            self._divgrad_matrix_sigmad = None
        if self.excitations.is_time_dependent:
            self.excitations.update(name="time", val=time)
            self._load_vector_electric = None

    @property
    def time(self) -> float:
        """The time instance of the solution, default = 0."""
        if self._time is None:
            return 0
        return self._time

    @time.setter
    def time(self, time: float):
        self._time = time
        self.ensure_consistency_time(time)

    @property
    def divgrad_matrix_eps(self) -> sparse.coo_matrix:
        """The divgrad matrix for the permittivity."""
        if self._divgrad_matrix_eps is None:
            self._divgrad_matrix_eps = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                            material.Permittivity)
        return self._divgrad_matrix_eps

    @divgrad_matrix_eps.setter
    def divgrad_matrix_eps(self, divgrad: sparse.coo_matrix):
        self._divgrad_matrix_eps = divgrad

    @property
    def divgrad_matrix_epsd(self) -> sparse.coo_matrix:
        """The divgrad matrix for the differential permittivity."""
        if self._divgrad_matrix_epsd is None:
            self._divgrad_matrix_epsd = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                             material.DifferentialPermittivity)
        return self._divgrad_matrix_epsd

    @divgrad_matrix_epsd.setter
    def divgrad_matrix_epsd(self, divgrad: sparse.coo_matrix):
        self._divgrad_matrix_epsd = divgrad

    @classmethod
    def extend_solution(cls, solution: 'CurrentFlowSolutionAxiStatic') -> 'CurrentFlowSolutionAxiIntermediate':
        """Creates an intermediate static solution from a basic static solution.

        Parameters
        ----------
        solution : CurrentFlowSolutionAxiStatic
            A basic stationary current solution.

        Returns
        -------
        intermediate_solution : CurrentFlowSolutionAxiIntermediate
        """
        intermediate_solution = cls.__call__(solution.description, solution.potential, solution.mesh,
                                             solution.shape_function, solution.regions, solution.materials,
                                             solution.excitations)
        return intermediate_solution


class CurrentFlowSolutionAxiHarmonic(HarmonicSolution):
    """The solution of a harmonic electroquasistatic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemAxiHarmonic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class CurrentFlowSolutionAxiTransient(TransientSolution):
    """The solution of a transient electroquasistatic problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.CurrentFlowProblemAxiTransient`.
    """

    related_static_solution = CurrentFlowSolutionAxiIntermediate

    _property_attributes = ['_nodes_not_on_axis', '_mass_one']

    _derived_attributes = ['_e_field', '_current_density',
                           '_joule_loss_power', '_joule_loss_power_density',
                           '_joule_loss_power_in_regions']

    _matrix_attributes = ['_divgrad_matrix_sigma', '_divgrad_matrix_sigmad',
                          '_divgrad_matrix_eps', '_divgrad_matrix_epsd',
                          '_load_vector_eqs']

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, potential: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_nodal_shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction',
                 regions: region.Regions,
                 materials: material.Materials, excitations: excitation.Excitations, time_steps: np.ndarray,
                 monitor_solutions: Dict[str, Tuple[np.ndarray, np.ndarray]]):
        """Constructor.

        Parameters
        ----------
        description :
        vector_potential :
        axi_mesh :
        tri_axi_nodal_shape_function :
        regions :
        materials :
        excitations :
        """
        super().__init__(description, potential, axi_mesh, tri_axi_nodal_shape_function, regions, materials,
                         excitations, time_steps, monitor_solutions)
        self.mesh: 'mesh.AxiMesh' = axi_mesh
        self.shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction' = tri_axi_nodal_shape_function
        self.plotter = TransientFieldPlotter(axi_mesh)

        self.consistency_check()

    def consistency_check(self):
        if not self.potential().shape == (len(self.time_steps), self.mesh.num_node):
            if self.potential().shape == (self.mesh.num_node, len(self.time_steps)):
                logger.warning("You passed the potential in wrong format. It is going to be transformed.")
                self.set_potential(self.potential().transpose())
            else:
                logger.critical("The potential has a unexpected format. It can't be used to compute other data.")
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

    def potential(self, *, time=None, index=None) -> np.ndarray:
        """The scalar potential. These are the coefficients for the basis functions.

        Parameters
        ----------
        time : float
            A time instance.
        index : int
            The index of a time instance.

        Returns
        -------
        potential : np.ndarray
        """
        return self._get_solution(time=time, index=index)

    def set_potential(self, potential: np.ndarray, *, time=None, index=None):
        """Set the potential.

        Parameters
        ----------
        potential : np.ndarray
            The potential to set
        time : float
            A time instance.
        index : int
            The index of a time instance.
        """
        self._set_solution(potential, time=time, index=index)

    # region Derived values
    @get_field('_e_field', 'electric field', 'ndarray')
    def e_field(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        return - self.shape_function.gradient(self.potential(index=index))

    @set_field('_e_field', 'electric field', 'ndarray')
    def set_e_field(self, e_field: Union[np.ndarray, Dict[int, np.ndarray]], *,
                    time=None, index=None):
        self._check_set_types(e_field, np.ndarray)

    @get_field('_current_density', 'current density', 'ndarray')
    def current_density(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        current_density_per_element = self.mesh.value_per_element(self.regions, self.materials, material.Conductivity)
        return np.tile(np.atleast_2d(current_density_per_element).T, 2) * self.e_field(index=index)

    @set_field('_current_density', 'current density', 'ndarray')
    def set_current_density(self, current_density: Union[np.ndarray, Dict[int, np.ndarray]], *,
                            time=None, index=None):
        self._check_set_types(current_density, np.ndarray)

    @get_field('_joule_loss_power', 'loss power', 'float')
    def joule_loss_power(self, *, time=None, index=None) -> Union[float, Dict[int, float]]:
        potential = self.potential(index=index)
        divgrad_matrix_sigma = self.divgrad_matrix_sigma(index=index)
        return potential @ divgrad_matrix_sigma @ potential

    @set_field('_joule_loss_power', 'loss power', 'float')
    def set_joule_loss_power(self, joule_loss_power: Union[float, Dict[int, float]], *, time=None, index=None):
        self._check_set_types(joule_loss_power, float)

    @get_field('_joule_loss_power_density', 'loss power density', 'ndarray')
    def joule_loss_power_density(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        e_field = self.e_field(index=index)
        current_density = self.current_density(index=index)
        return np.sum(e_field * current_density, axis=1)

    @set_field('_joule_loss_power_density', 'loss power density', 'ndarray')
    def set_joule_loss_power_density(self, joule_loss_power_density: Union[np.ndarray, Dict[int, np.ndarray]], *,
                                     time=None, index=None):
        self._check_set_types(joule_loss_power_density, np.ndarray)

    def joule_loss_power_in_regions(self, *region_ids: int, time=None, index=None) -> float:
        """Compute the loss power in certain regions.

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
        # pylint: disable=consider-using-generator
        if time is None and index is None:
            return self._joule_loss_power_in_regions
        index = self.get_index(time, index)

        # If the time index has not calculated yet, add an empty dict
        if index not in self._joule_loss_power_in_regions.keys():
            self._joule_loss_power_in_regions[index] = {}

        indices_per_region = []  # List of indices of the elements in the regions
        for region_id in region_ids:
            if not self.regions.get_regi(region_id).dim == 2:
                raise ValueError(f"The region '{region_id}' does not belong to a region of the expected dimension.")
            if region_id not in self._joule_loss_power_in_regions[index].keys():
                indices = indices_per_region.append(np.where(self.mesh.elem2regi == region_id)[0])
                self._joule_loss_power_in_regions[index][region_id] = 2 * np.pi * np.sum(
                    self.joule_loss_power_density(index=index)[indices] * self.mesh.elem_area[indices] *
                    self.mesh.node[0, :])
        return sum([self._joule_loss_power_in_regions[index][region_id] for region_id in region_ids])

    # endregion

    # region Matrix and vector

    @get_matrix('_divgrad_matrix_sigma', 'divgrad_matrix (Conductivity)')
    def divgrad_matrix_sigma(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.Conductivity)
        time_dependent = self.materials.is_time_dependent(material.Conductivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.Conductivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.Conductivity)

        return lin and not time_dependent, lambda: self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                                        material.Conductivity)

    @set_matrix('_divgrad_matrix_sigma', 'divgrad matrix (Conductivity)')
    def set_divgrad_matrix_sigma(self, divgrad_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_divgrad_matrix_sigmad', 'divgrad_matrix (DifferentialConductivity)')
    def divgrad_matrix_sigmad(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.DifferentialConductivity)
        time_dependent = self.materials.is_time_dependent(material.DifferentialConductivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.DifferentialConductivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.DifferentialConductivity)

        return lin and not time_dependent, \
               lambda: self.shape_function.divgrad_operator(self.regions, self.materials,
                                                            material.DifferentialConductivity)

    @set_matrix('_divgrad_matrix_sigmad', 'divgrad matrix (DifferentialConductivity)')
    def set_divgrad_matrix_sigmad(self, divgrad_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_divgrad_matrix_eps', 'divgrad_matrix (Permittivity)')
    def divgrad_matrix_eps(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.Permittivity)
        time_dependent = self.materials.is_time_dependent(material.Permittivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.Permittivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.Permittivity)

        return lin and not time_dependent, lambda: self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                                        material.Permittivity)

    @set_matrix('_divgrad_matrix_eps', 'divgrad matrix (Permittivity)')
    def set_divgrad_matrix_eps(self, divgrad_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_divgrad_matrix_epsd', 'divgrad_matrix (DifferentialPermittivity)')
    def divgrad_matrix_epsd(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.DifferentialPermittivity)
        time_dependent = self.materials.is_time_dependent(material.DifferentialPermittivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.DifferentialPermittivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.DifferentialPermittivity)

        return lin and not time_dependent, \
               lambda: self.shape_function.divgrad_operator(self.regions,
                                                            self.materials,
                                                            material.DifferentialPermittivity)

    @set_matrix('_divgrad_matrix_epsd', 'divgrad matrix (DifferentialPermittivity)')
    def set_divgrad_matrix_epsd(self, divgrad_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_load_vector_eqs', 'load vector (EQS)')
    def load_vector_eqs(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.excitations.is_linear
        time_dependent = self.excitations.is_time_dependent

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.excitations.update('solution', static_solution)
        if time_dependent:
            self.excitations.update('time', self.time_steps[index])

        return lin and not time_dependent, lambda: self.shape_function.load_vector(self.regions, self.excitations)

    @set_matrix('_load_vector_eqs', 'load vector (EQS)')
    def set_load_vector_eqs(self, load_vector: sparse.coo_matrix, *, time=None, index=None):
        pass

    # endregion

    # region Plot on time instance and animate

    def plot_potential(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the potential.

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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'potential in V'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.potential(index=index),
                                              self.time_steps[index], **kwargs)

    def animate_potential(self, *, times: Iterable[float] = None,
                          indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the potential.

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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'potential in V'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.animate_scalar_field((self.potential(index=i) for i in indices),
                                                 **kwargs)

    def plot_e_field(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> \
            Tuple[Figure, Type[Axes]]:
        """Plot the electric field.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Electric field in V/m"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_vector_field(self.e_field(index=index),
                                              plot_type, self.time_steps[index], **kwargs)

    def animate_e_field(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                        indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the electric field.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Electric field in V/m"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.e_field(index=i) for i in indices), plot_type, **kwargs)

    def plot_current_density(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> \
            Tuple[Figure, Type[Axes]]:
        """Plot the current density.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Current density in A/m^2"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.current_density(index=index), plot_type, self.time_steps[index],
                                              **kwargs)

    def animate_current_density(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                                indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the current density.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Current density in A/m^2"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.current_density(index=i) for i in indices), plot_type, **kwargs)

    def plot_joule_loss_power_density(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis',
                          'title': 'Loss power density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.joule_loss_power_density(index=index),
                                              self.time_steps[index], **kwargs)

    def animate_joule_loss_power_density(self, *, times: Iterable[float] = None,
                                         indices: Union[Iterable[int], int] = None,
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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"loss power density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.joule_loss_power_density(index=i) for i in indices), **kwargs)

    def plot_load_vector_eqs(self, *, time=None, index=None, **kwargs):
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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'Load vector'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.load_vector_eqs(index=index).toarray(), self.time_steps[index],
                                              **kwargs)

    def animate_load_vector_eqs(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
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
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Load vector"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.load_vector_eqs(index=i).toarray() for i in indices), **kwargs)

    def plot_equilines_potential(self, *, time=None, index=None, **kwargs):
        """Plot the equilines of the potential.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'Equilines potential',
                          'colors': 'b'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_equilines(self.potential(index=index),
                                           self.time_steps[index], **kwargs)

    def animate_equilines_potential(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                                    **kwargs):
        """Animate the equilines of the potential.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': "Equilines potential"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.potential(index=i) for i in indices), **kwargs)

    # endregion

    def plot_monitor_solution(self, key: str, **kwargs):
        if key not in self.monitor_solutions:
            raise KeyError(f"Key '{key}' not known in monitor solutions")
        return self.plotter.plot_over_time(self.monitor_solutions[key][0], self.monitor_solutions[key][1], **kwargs)
