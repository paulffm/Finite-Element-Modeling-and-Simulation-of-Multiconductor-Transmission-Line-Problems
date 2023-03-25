# coding=utf-8
"""Thermal solutions in axisymmetric coordinates

.. sectionauthor:: Ruppert, Bundschuh
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

from pyrit.solution.Solution import StaticSolution, HarmonicSolution, TransientSolution, get_field, set_field, \
    get_matrix, set_matrix
from pyrit.solution.FieldPlotter import StaticFieldPlotter, TransientFieldPlotter

if TYPE_CHECKING:
    from pyrit import shapefunction, mesh

logger = get_logger(__name__)

__all__ = ['ThermalSolutionAxiStatic', 'ThermalSolutionAxiHarmonic', 'ThermalSolutionAxiTransient']


class ThermalSolutionAxiStatic(StaticSolution):
    """The solution of a stationary heat conduction problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemAxiStatic`.
    """

    #: Protected attribute names that are used for properties
    _derived_attributes = ['_temperature_gradient', '_temperature_gradient_abs',
                           '_heat_flux_density', '_heat_flux_density_abs',
                           '_thermal_conductivity_per_element', '_thermal_loss_power',
                           '_thermal_loss_power_density', '_thermal_loss_power_in_regions']

    _matrix_attributes = ['_divgrad_matrix_lambda', '_load_vector_thermal']

    _property_attributes = _derived_attributes.copy() + _matrix_attributes.copy()

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, temperature: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_nodal_shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction',
                 regions: region.Regions, materials: material.Materials, excitations: excitation.Excitations):
        """The solution of a stationary heat conduction problem in axisymmetric coordinates.

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
        super().__init__(description, temperature, None, None, regions, materials, excitations)
        self.mesh: mesh.AxiMesh = axi_mesh
        self.shape_function: shapefunction.TriAxisymmetricNodalShapeFunction = tri_axi_nodal_shape_function
        self.plotter = StaticFieldPlotter(axi_mesh)

        self.consistency_check()
        self.update_material_thermal()

    def update_material_thermal(self):
        # Write solution to nonlinear materials and excitations to ensure consistency.
        self.materials.update(name="solution", val=self, prop_class=material.ThermalConductivity)
        self.materials.update(name="solution", val=self, prop_class=material.VolumetricHeatCapacity)
        self.materials.update(name="temperature", val=self.temperature, prop_class=material.Conductivity)
        self.materials.update(name="temperature", val=self.temperature, prop_class=material.DifferentialConductivity)
        self.materials.update(name="temperature", val=self.temperature, prop_class=material.Permittivity)
        self.materials.update(name="temperature", val=self.temperature, prop_class=material.DifferentialPermittivity)
        self.excitations.update(name="solution", val=self)

    def ensure_consistency_temperature(self):
        # Function that ensures consistency by resetting all values derived from the temperature.
        for prop_name in ThermalSolutionAxiStatic._derived_attributes:
            setattr(self, prop_name, None)

        self.update_material_thermal()

        if not self.materials.is_linear(material.ThermalConductivity):
            self._divgrad_matrix_lambda = None
        if not self.excitations.is_linear:
            self._load_vector_thermal = None

    @property
    def temperature(self) -> np.ndarray:
        """The temperature. These are the coefficients for the basis functions."""
        return self._solution

    @temperature.setter
    def temperature(self, temperature: np.ndarray):
        self._solution = temperature
        self.ensure_consistency_temperature()

    @property
    def temperature_gradient(self) -> np.ndarray:
        """The temperature gradient. Array with two values per mesh element."""
        if self._temperature_gradient is None:
            self._temperature_gradient = - self.shape_function.gradient(self.temperature)
        return self._temperature_gradient

    @temperature_gradient.setter
    def temperature_gradient(self, temperature_gradient: np.ndarray):
        self._temperature_gradient = temperature_gradient

    @property
    def temperature_gradient_abs(self) -> np.ndarray:
        """Absolute value of the temperature gradient. Array with one value per mesh element."""
        if self._temperature_gradient_abs is None:
            self._temperature_gradient_abs = np.linalg.norm(self.temperature_gradient, axis=1)
        return self._temperature_gradient_abs

    @temperature_gradient_abs.setter
    def temperature_gradient_abs(self, temperature_gradient_abs: np.ndarray):
        self._temperature_gradient_abs = temperature_gradient_abs

    @property
    def thermal_conductivity_per_element(self) -> np.ndarray:
        """The thermal conductivity. Array with one value per mesh element."""
        if self._thermal_conductivity_per_element is None:
            self._thermal_conductivity_per_element = self.mesh.value_per_element(self.regions,
                                                                                 self.materials,
                                                                                 material.ThermalConductivity)
        return self._thermal_conductivity_per_element

    @thermal_conductivity_per_element.setter
    def thermal_conductivity_per_element(self, thermal_cond: np.ndarray):
        self._thermal_conductivity_per_element = thermal_cond

    @property
    def heat_flux_density(self) -> np.ndarray:
        """The heat flux density. Array with two values per mesh element."""
        if self._heat_flux_density is None:
            self._heat_flux_density = self.thermal_conductivity_per_element[:, None] * self.temperature_gradient
        return self._heat_flux_density

    @heat_flux_density.setter
    def heat_flux_density(self, heat_flux_density):
        self._heat_flux_density = heat_flux_density

    @property
    def heat_flux_density_abs(self) -> np.ndarray:
        """The absolute value of the heat flux density. Array with two values per mesh element."""
        if self._heat_flux_density_abs is None:
            self._heat_flux_density_abs = self.thermal_conductivity_per_element * self.temperature_gradient_abs
        return self._heat_flux_density_abs

    @heat_flux_density_abs.setter
    def heat_flux_density_abs(self, heat_flux_density_abs):
        self._heat_flux_density_abs = heat_flux_density_abs

    @property
    def divgrad_matrix_lambda(self) -> sparse.coo_matrix:
        """The divgrad matrix for the thermal conductivity."""
        if self._divgrad_matrix_lambda is None:
            self._divgrad_matrix_lambda = self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                               material.ThermalConductivity)
        return self._divgrad_matrix_lambda

    @divgrad_matrix_lambda.setter
    def divgrad_matrix_lambda(self, divgrad: sparse.coo_matrix):
        self._divgrad_matrix_lambda = divgrad

    @property
    def load_vector_thermal(self) -> sparse.coo_matrix:
        """The load vector of the problem."""
        if self._load_vector_thermal is None:
            self._load_vector_thermal = self.shape_function.load_vector(self.regions, self.excitations)
        return self._load_vector_thermal

    @load_vector_thermal.setter
    def load_vector_thermal(self, load: sparse.coo_matrix):
        self._load_vector_thermal = load

    @property
    def thermal_loss_power_density(self) -> np.ndarray:
        """The thermal loss power density. Array with one value per mesh element."""
        if self._thermal_loss_power_density is None:
            self._thermal_loss_power_density = np.sum(self.heat_flux_density * self.temperature_gradient, axis=1)
        return self._thermal_loss_power_density

    @thermal_loss_power_density.setter
    def thermal_loss_power_density(self, pld):
        self._thermal_loss_power_density = pld

    # ToDo: Check formular for thermal_loss_power and change name. Is there a 0.5 missing?
    @property
    def thermal_loss_power(self) -> float:
        """The thermal loss power."""
        if self._thermal_loss_power is None:
            self._thermal_loss_power = self.temperature @ self.divgrad_matrix_lambda @ self.temperature
        return self._thermal_loss_power

    @thermal_loss_power.setter
    def thermal_loss_power(self, thermal_loss_power: float):
        self._thermal_loss_power = thermal_loss_power

    def thermal_loss_power_in_regions(self, *region_ids: int) -> float:
        """Compute the thermal power loss in certain regions.

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
        energy = 2 * np.pi * np.sum(self.thermal_loss_power_density[indices] *
                                    self.mesh.elem_area[indices] * mean_radii)

        return energy

    def consistency_check(self):
        if self.temperature.shape != (self.mesh.num_node,):
            logger.critical("The temperature has an unexpected format. It can't be used to compute other data.")

    def plot_temperature(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the temperature.

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
                          'title': 'temperature in K'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.temperature, **kwargs)

    def plot_temperature_gradient(self, kind: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the temperature gradient.

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
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"temperature gradient"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        if kind == 'line':
            kwargs.setdefault('colors', 'blue')
            kwargs.setdefault('num_lines', 20)
            return self.plotter.plot_equilines(self.temperature, **kwargs)
        return self.plotter.plot_vector_field(self.temperature_gradient, kind, **kwargs)

    def plot_heat_flux_density(self, kind: str = 'arrows', **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the heat flux density.

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
        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"heat flux density in W/m^2"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.heat_flux_density, kind, **kwargs)

    def plot_thermal_loss_power_density(self, **kwargs) -> Tuple[Figure, Type[Axes]]:
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
        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis',
                          'title': 'thermal power loss density W/m^3'}
        default_kwargs.update(kwargs)
        return self.plotter.plot_scalar_field(self.thermal_loss_power_density, **default_kwargs)


class ThermalSolutionAxiIntermediate(ThermalSolutionAxiStatic):
    """An intermediate solution of a thermal problem in axisymmetric coordinates.

    An intermediate solution of a thermal problem. It has all properties of a stationary heat conduction solution, but
    is extended by properties concerning time-dependencies and the volumetric heat capacity. It is used for internal
    auxiliary variables in the solve-method of :py:class:`~pyrit.problem.problems.ThermalProblemAxiTransient`.
    """

    #: Protected attribute names that are used for properties
    _derived_attributes = ['_temperature_gradient', '_temperature_gradient_abs',
                           '_heat_flux_density', '_heat_flux_density_abs',
                           '_thermal_conductivity_per_element', '_thermal_loss_power',
                           '_thermal_loss_power_density', '_thermal_loss_power_in_regions']

    _matrix_attributes = ['_divgrad_matrix_lambda', '_mass_matrix_s', '_load_vector_thermal']

    _property_attributes = ['_time'] + _derived_attributes.copy() + _matrix_attributes.copy()

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, temperature: np.ndarray, axi_mesh: 'mesh.AxiMesh',
                 tri_axi_nodal_shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction',
                 regions: region.Regions, materials: material.Materials, excitations: excitation.Excitations):
        """The solution of a stationary heat conduction problem in axisymmetric coordinates.

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
        super().__init__(description, temperature, axi_mesh, tri_axi_nodal_shape_function, regions,
                         materials, excitations)

    def ensure_consistency_temperature(self):
        # Function that ensures consistency by resetting all values derived from the temperature.
        super().ensure_consistency_temperature()

        if not self.materials.is_linear(material.VolumetricHeatCapacity):
            self._mass_matrix_s = None

    def ensure_consistency_time(self, time):
        # Write time to time-dependent materials and excitations to ensure consistency.
        if self.materials.is_time_dependent(prop_class=material.ThermalConductivity):
            self.materials.update(name="time", val=time, prop_class=material.ThermalConductivity)
            self._divgrad_matrix_lambda = None
        if self.materials.is_time_dependent(prop_class=material.VolumetricHeatCapacity):
            self.materials.update(name="time", val=time, prop_class=material.VolumetricHeatCapacity)
            self._mass_matrix_s = None
        if self.excitations.is_time_dependent:
            self.excitations.update(name="time", val=time)
            self._load_vector_thermal = None

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
    def mass_matrix_s(self) -> sparse.coo_matrix:
        """The mass matrix for the volumetric heat capacity."""
        if self._mass_matrix_s is None:
            self._mass_matrix_s = self.shape_function.mass_matrix(self.regions, self.materials,
                                                                  material.VolumetricHeatCapacity)
        return self._mass_matrix_s

    @mass_matrix_s.setter
    def mass_matrix_s(self, mass: sparse.coo_matrix):
        self._mass_matrix_s = mass

    @classmethod
    def extend_solution(cls, solution: 'ThermalSolutionAxiStatic') -> 'ThermalSolutionAxiIntermediate':
        """Creates an intermediate static solution from a basic static solution.

        Parameters
        ----------
        solution : ThermalSolutionAxiStatic
            A basic static thermal solution.

        Returns
        -------
        intermediate_solution : ThermalSolutionAxiIntermediate
        """
        intermediate_solution = cls.__call__(solution.description, solution.temperature, solution.mesh,
                                             solution.shape_function, solution.regions, solution.materials,
                                             solution.excitations)
        return intermediate_solution


class ThermalSolutionAxiHarmonic(HarmonicSolution):
    """The solution of a harmonic heat conduction problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemAxiHarmonic`.
    """

    def consistency_check(self):
        raise NotImplementedError()


class ThermalSolutionAxiTransient(TransientSolution):
    """The solution of a transient heat conduction problem in axisymmetric coordinates.

    The corresponding problem class is :py:class:`~pyrit.problem.problems.ThermalProblemAxiTransient`.
    """

    related_static_solution = ThermalSolutionAxiIntermediate

    _property_attributes = ['_nodes_not_on_axis']

    _derived_attributes = ['_temperature_gradient', '_heat_flux_density',
                           '_thermal_loss_power', '_thermal_loss_power_density',
                           '_thermal_loss_power_in_regions']

    _matrix_attributes = ['_divgrad_matrix_lambda', '_mass_matrix_s', '_load_vector_thermal']

    ignore_for_saving = _property_attributes.copy() + _derived_attributes.copy() + _matrix_attributes.copy()

    def __init__(self, description: str, temperature: np.ndarray, axi_mesh: 'mesh.AxiMesh',
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
        super().__init__(description, temperature, axi_mesh, tri_axi_nodal_shape_function, regions, materials,
                         excitations, time_steps, monitor_solutions)
        self.mesh: 'mesh.AxiMesh' = axi_mesh
        self.shape_function: 'shapefunction.TriAxisymmetricNodalShapeFunction' = tri_axi_nodal_shape_function
        self.plotter = TransientFieldPlotter(axi_mesh)

        self.consistency_check()

    def consistency_check(self):
        if not self.temperature().shape == (len(self.time_steps), self.mesh.num_node):
            if self.temperature().shape == (self.mesh.num_node, len(self.time_steps)):
                logger.warning("You passed the temperature in wrong format. It is going to be transformed.")
                self.set_temperature(self.temperature().transpose())
            else:
                logger.critical("The temperature has a unexpected format. It can't be used to compute other data.")
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

    def temperature(self, *, time=None, index=None) -> np.ndarray:
        """The temperature. These are the coefficients for the basis functions.

        Parameters
        ----------
        time : float
            A time instance.
        index : int
            The index of a time instance.

        Returns
        -------
        temperature : np.ndarray
        """
        return self._get_solution(time=time, index=index)

    def set_temperature(self, temperature: np.ndarray, *, time=None, index=None):
        """Set the temperature.

        Parameters
        ----------
        temperature : np.ndarray
            The temperature to set
        time : float
            A time instance.
        index : int
            The index of a time instance.
        """
        self._set_solution(temperature, time=time, index=index)

    # region Derived values
    @get_field('_temperature_gradient', 'temperature gradient', 'ndarray')
    def temperature_gradient(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        return - self.shape_function.gradient(self.temperature(index=index))

    @set_field('_temperature_gradient', 'temperature gradient', 'ndarray')
    def set_temperature_gradient(self, temperature_gradient: Union[np.ndarray, Dict[int, np.ndarray]], *,
                                 time=None, index=None):
        self._check_set_types(temperature_gradient, np.ndarray)

    @get_field('_heat_flux_density', 'heat flux density', 'ndarray')
    def heat_flux_density(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        thermal_conductivity_per_element = self.mesh.value_per_element(self.regions,
                                                                       self.materials, material.ThermalConductivity)
        return np.tile(np.atleast_2d(thermal_conductivity_per_element).T, 2) * self.temperature_gradient(index=index)

    @set_field('_heat_flux_density', 'heat flux density', 'ndarray')
    def set_heat_flux_density(self, heat_flux_density: Union[np.ndarray, Dict[int, np.ndarray]], *,
                              time=None, index=None):
        self._check_set_types(heat_flux_density, np.ndarray)

    @get_field('_thermal_loss_power', 'thermal loss power', 'float')
    def thermal_loss_power(self, *, time=None, index=None) -> Union[float, Dict[int, float]]:
        temperature = self.temperature(index=index)
        divgrad_matrix = self.divgrad_matrix_lambda(index=index)
        return temperature @ divgrad_matrix @ temperature

    @set_field('_thermal_loss_power', 'thermal loss power', 'float')
    def set_thermal_loss_power(self, thermal_loss_power: Union[float, Dict[int, float]], *, time=None, index=None):
        self._check_set_types(thermal_loss_power, float)

    @get_field('_thermal_loss_power_density', 'thermal loss power density', 'ndarray')
    def thermal_loss_power_density(self, *, time=None, index=None) -> Union[np.ndarray, Dict[int, np.ndarray]]:
        temperature_gradient = self.temperature_gradient(index=index)
        heat_flux_density = self.heat_flux_density(index=index)
        return np.sum(temperature_gradient * heat_flux_density, axis=1)

    @set_field('_thermal_loss_power_density', 'thermal loss power density', 'ndarray')
    def set_thermal_loss_power_density(self, thermal_loss_power_density: Union[np.ndarray, Dict[int, np.ndarray]], *,
                                       time=None, index=None):
        self._check_set_types(thermal_loss_power_density, np.ndarray)

    def thermal_loss_power_in_regions(self, *region_ids: int, time=None, index=None) -> float:
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
        # pylint: disable=consider-using-generator
        if time is None and index is None:
            return self._energy_in_regions
        index = self.get_index(time, index)

        # If the time index has not calculated yet, add an empty dict
        if index not in self._energy_in_regions.keys():
            self._thermal_loss_power_in_regions[index] = {}

        indices_per_region = []  # List of indices of the elements in the regions
        for region_id in region_ids:
            if not self.regions.get_regi(region_id).dim == 2:
                raise ValueError(f"The region '{region_id}' does not belong to a region of the expected dimension.")
            if region_id not in self._thermal_loss_power_in_regions[index].keys():
                indices = indices_per_region.append(np.where(self.mesh.elem2regi == region_id)[0])
                self._thermal_loss_power_in_regions[index][region_id] = 2 * np.pi * np.sum(
                    self.thermal_loss_power_density(index=index)[indices] * self.mesh.elem_area[indices] *
                    self.mesh.node[0, :])

        return sum([self._thermal_loss_power_in_regions[index][region_id] for region_id in region_ids])

    # endregion

    # region Matrix and vector

    @get_matrix('_divgrad_matrix_lambda', 'divgrad matrix (thermal)')
    def divgrad_matrix_lambda(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.ThermalConductivity)
        time_dependent = self.materials.is_time_dependent(material.ThermalConductivity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.ThermalConductivity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.ThermalConductivity)

        return lin and not time_dependent, lambda: self.shape_function.divgrad_operator(self.regions, self.materials,
                                                                                        material.ThermalConductivity)

    @set_matrix('_divgrad_matrix_lambda', 'divgrad matrix (thermal)')
    def set_divgrad_matrix_lambda(self, divgrad_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_mass_matrix_s', 'mass matrix (thermal)')
    def mass_matrix_s(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.materials.is_linear(material.VolumetricHeatCapacity)
        time_dependent = self.materials.is_time_dependent(material.VolumetricHeatCapacity)

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.materials.update('solution', static_solution, prop_class=material.VolumetricHeatCapacity)
        if time_dependent:
            self.materials.update('time', self.time_steps[index], prop_class=material.VolumetricHeatCapacity)

        return lin and not time_dependent, lambda: self.shape_function.mass_matrix(self.regions, self.materials,
                                                                                   material.VolumetricHeatCapacity)

    @set_matrix('_mass_matrix_s', 'mass matrix (thermal)')
    def set_mass_matrix_s(self, mass_matrix: sparse.coo_matrix, *, time=None, index=None):
        pass

    @get_matrix('_load_vector_thermal', 'load vector (thermal)')
    def load_vector_thermal(self, *, time=None, index=None) -> Union[Dict[int, sparse.coo_matrix], sparse.coo_matrix]:
        lin = self.excitations.is_linear
        time_dependent = self.excitations.is_time_dependent

        if not lin:
            static_solution = self.get_static_solution(index=index)
            self.excitations.update('solution', static_solution)
        if time_dependent:
            self.excitations.update('time', self.time_steps[index])

        return lin and not time_dependent, lambda: self.shape_function.load_vector(self.regions, self.excitations)

    @set_matrix('_load_vector_thermal', 'load vector (thermal)')
    def set_load_vector_thermal(self, load_vector: sparse.coo_matrix, *, time=None, index=None):
        pass

    # endregion

    # region Plot on time instance and animate

    def plot_temperature(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
        """Plot the temperature.

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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'temperature in K'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.temperature(index=index),
                                              self.time_steps[index], **kwargs)

    def animate_temperature(self, *, times: Iterable[float] = None,
                            indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the temperature.

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

        default_kwargs = {'plot_3d': False, 'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'temperature in K'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.animate_scalar_field((self.temperature(index=i) for i in indices),
                                                 **kwargs)

    def plot_heat_flux_density(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> \
            Tuple[Figure, Type[Axes]]:
        """Plot the heat flux density.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Heat flux density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_vector_field(self.heat_flux_density(index=index),
                                              plot_type, self.time_steps[index], **kwargs)

    def animate_heat_flux_density(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                                  indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the heat flux density.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"Heat flux density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.heat_flux_density(index=i) for i in indices), plot_type, **kwargs)

    def plot_temperature_gradient(self, plot_type: str = 'arrows', *, time=None, index=None, **kwargs) -> \
            Tuple[Figure, Type[Axes]]:
        """Plot the temperature gradient.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"temperature gradient"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

        return self.plotter.plot_vector_field(self.temperature_gradient(index=index), plot_type, self.time_steps[index],
                                              **kwargs)

    def animate_temperature_gradient(self, plot_type: str = 'arrows', *, times: Iterable[float] = None,
                                     indices: Union[Iterable[int], int] = None, **kwargs):
        """Animate the temperature gradient.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"temperature gradient"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_vector_field((self.temperature_gradient(index=i) for i in indices), plot_type, **kwargs)

    def plot_thermal_loss_power_density(self, *, time=None, index=None, **kwargs) -> Tuple[Figure, Type[Axes]]:
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
                          'title': 'Thermal loss power density'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_scalar_field(self.thermal_loss_power_density(index=index),
                                              self.time_steps[index], **kwargs)

    def animate_thermal_loss_power_density(self, *, times: Iterable[float] = None,
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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': r"thermal loss power density"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.thermal_loss_power_density(index=i) for i in indices), **kwargs)

    def plot_load_vector_thermal(self, *, time=None, index=None, **kwargs):
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
        return self.plotter.plot_scalar_field(self.load_vector_thermal(index=index).toarray(), self.time_steps[index],
                                              **kwargs)

    def animate_load_vector_thermal(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
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
        self.plotter.animate_scalar_field((self.load_vector_thermal(index=i).toarray() for i in indices), **kwargs)

    def plot_equilines_temperature(self, *, time=None, index=None, **kwargs):
        """Plot the equilines of the temperature.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': 'Equilines temperature',
                          'colors': 'b'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        return self.plotter.plot_equilines(self.temperature(index=index),
                                           self.time_steps[index], **kwargs)

    def animate_equilines_temperature(self, *, times: Iterable[float] = None, indices: Union[Iterable[int], int] = None,
                                      **kwargs):
        """Animate the equilines of the temperature.

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

        default_kwargs = {'x_label': 'rho axis', 'y_label': 'z axis', 'title': "Equilines temperature"}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)
        self.plotter.animate_scalar_field((self.temperature(index=i) for i in indices), **kwargs)

    # endregion
    def plot_monitor_solution(self, key: str, **kwargs):
        if key not in self.monitor_solutions:
            raise KeyError(f"Key '{key}' not known in monitor solutions")
        return self.plotter.plot_over_time(self.monitor_solutions[key][0], self.monitor_solutions[key][1], **kwargs)
