# coding=utf-8
"""Current flow problem in axisymmetric coordinates

.. sectionauthor:: Ruppert, Bundschuh
"""
# pylint: disable=arguments-differ
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, Union, Tuple, Callable, NoReturn, Any, Iterable
from scipy.sparse import coo_matrix
import numpy as np

from pyrit import get_logger
from pyrit import mesh, shapefunction, region, material, bdrycond, excitation
from pyrit.excitation import Excitations
from pyrit.problem.Problem import StaticProblem, HarmonicProblem, TransientProblem
from pyrit.solution import CurrentFlowSolutionAxiStatic, CurrentFlowSolutionAxiHarmonic, \
    CurrentFlowSolutionAxiTransient
from pyrit.solution.CurrentFlowSolutionAxi import CurrentFlowSolutionAxiIntermediate

if TYPE_CHECKING:
    from pyrit.mesh import AxiMesh
    from pyrit.shapefunction import TriAxisymmetricNodalShapeFunction

logger = get_logger(__name__)

__all__ = ['CurrentFlowProblemAxiStatic', 'CurrentFlowProblemAxiHarmonic', 'CurrentFlowProblemAxiTransient',
           'CurrentFlowSolverInfoAxiStatic', 'CurrentFlowSolverInfoAxiTransient']


@dataclass()  # todo this should be explained better/removed from the doc, as then it's for internal/expert use only.
class CurrentFlowSolverInfoAxiStatic:
    """Class for passing options to the solve-method of CurrentFlowProblemAxiStatic"""

    tolerance_newton: float = 1e-8  # tolerance for the newton algorithm
    max_iter_newton: int = 100  # maximum number of newton iterations
    relaxation_newton: float = 1  # relaxation for newton algorithm


@dataclass()
class CurrentFlowSolverInfoAxiTransient:
    """Class for passing options to the solve-method of CurrentFlowProblemAxiTransient"""

    tolerance_newton: float = 1e-8  # tolerance for the newton algorithm
    max_iter_newton: int = 100  # maximum number of newton iterations
    relaxation_newton: float = 1  # relaxation for newton algorithm


class CurrentFlowProblemAxiStatic(StaticProblem):
    r"""A stationary current problem in axisymmetric coordinates:

    The stationary current problem models resistive effects. The corresponding differential equation reads

    .. math::
        -\mathrm{div}(\sigma \, \mathrm{grad} (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`~pyrit.material.Conductivity`) and :math:`\phi`
    denotes the electric scalar potential. A possible application is, for example, the steady state simulation of a HVDC
    cable joint.

    In case of nonlinear problems, the solve routine performs a Newton algorithm. The algorithm requires the definition
    of a :py:class:`~pyrit.material.DifferentialConductivity` for all materials with a nonlinear conductivity.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.CurrentFlowSolutionAxiStatic`.
    """

    problem_identifier: str = 'Stationary current problem in axisymmetric coordinates'

    def __init__(self, description: str, axi_mesh: mesh.AxiMesh,
                 tri_nodal_shape_function: shapefunction.TriAxisymmetricNodalShapeFunction,
                 regions: region.Regions, materials: material.Materials,
                 boundary_conditions: bdrycond.BdryCond,
                 excitations: excitation.Excitations = None):
        """Constructor

        Parameters
        ----------
        description : str
            A description of the problem
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`
        tri_nodal_shape_function : TriAxisymmetricNodalShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricNodalShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`
        boundary_conditions : BdryCond
            A boundary conditions object representing the boundary conditions for the stationary current problem.
            See :py:mod:`pyrit.bdrycond`
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`
        """
        # Setup Model
        super().__init__(description, None, None, regions, materials, boundary_conditions, excitations)
        # So the type is known by the IDE
        self.mesh: mesh.AxiMesh = axi_mesh
        self.shape_function: shapefunction.TriAxisymmetricNodalShapeFunction = tri_nodal_shape_function
        if excitations is None:
            self.excitations = Excitations()

        # Determine if the problem is linear (avoid Newton algorithm)
        self.is_linear: bool = materials.is_linear(prop_class=material.Conductivity) and \
                               self.excitations.is_linear and boundary_conditions.is_linear

        # Determine if problem is temperature-dependent (needed for coupled problems)
        self._is_temperature_dependent: bool = None
        self.consistency_check()

    @property
    def is_temperature_dependent(self):
        """Checks if the problem is temperature dependent."""
        if self._is_temperature_dependent is None:
            self._is_temperature_dependent = False
            for mat in self.materials:
                if 'temperature' in mat.get_property(material.Conductivity).keyword_args:
                    self._is_temperature_dependent = True
            for exci in self.excitations:
                if 'temperature' in exci.keyword_args:
                    self._is_temperature_dependent = True
            for bc in self.boundary_conditions:
                if 'temperature' in bc.keyword_args:
                    self._is_temperature_dependent = True
        return self._is_temperature_dependent

    def consistency_check(self):
        # make sure all materials with a nonlinear conductivity also have a differential conductivity
        if not self.is_linear:
            for mat_id in self.materials.get_ids():
                mat = self.materials.get_material(mat_id)
                cond = mat.get_property(material.Conductivity)
                if cond is not None:
                    if not cond.is_linear:
                        if mat.get_property(material.DifferentialConductivity) is None:
                            raise ValueError("All nonlinear materials must have a DifferentialConductivity. "
                                             "The DifferentialConductivity of " + mat.name + " is missing.")

    def _solve_system(self, matrix, rhs, **kwargs) -> np.ndarray:
        # Help routine that shrinks, solves and inflates the system "matrix * x = rhs", returns the solution x.

        matrix_shrink, rhs_shrink, _, _, support_data = self.shape_function.shrink(matrix.tocoo(), coo_matrix(rhs),
                                                                                   self, 1)

        potential_shrink, _ = type(self).solve_linear_system(matrix_shrink.tocsr(), rhs_shrink.todense(),
                                                             **kwargs)
        return self.shape_function.inflate(potential_shrink, self, support_data)

    def solve(self, u0: Union[np.ndarray, CurrentFlowSolutionAxiStatic] = None,
              solver_info: CurrentFlowSolverInfoAxiStatic = None,
              **kwargs) -> CurrentFlowSolutionAxiStatic:  # pylint: disable=arguments-differ
        """Solve the Problem using the Newton method.

        Parameters
        ----------
        u0: np.ndarray
            [V] Start guess for Newton iteration. If None is given, an all-zero vector is used.
        solver_info: CurrentFlowSolverInfoAxiStatic
            Options for the Newton algorithm
        kwargs :
            Are all passed to the function :py:func:`solve_linear_system`.

        Returns
        -------
        solution : CurrentFlowSolutionAxiStatic
            The solution object
        """
        # region Linear Probelm
        if self.is_linear:  # linear problem
            static_solution = CurrentFlowSolutionAxiStatic.solution_from_problem(self.description,
                                                                                 self, np.zeros((self.mesh.num_node,)))
            static_solution.potential = self._solve_system(static_solution.divgrad_matrix_sigma,
                                                           static_solution.load_vector_electric, **kwargs)
            return static_solution
        # endregion

        # region Nonlinear Problem --> Newton algorithm
        # region Initializations
        iter_no = 0
        if solver_info is None:
            solver_info = CurrentFlowSolverInfoAxiStatic()

        # Set initial guess for solution
        if u0 is None:
            # Initialize materials with zero solution and solve boundary value problem for this material
            # working point
            static_solution = CurrentFlowSolutionAxiStatic.solution_from_problem(self.description,
                                                                                 self, np.zeros((self.mesh.num_node,)))
            self.boundary_conditions.update('solution', static_solution)
            static_solution.potential = self._solve_system(static_solution.divgrad_matrix_sigma,
                                                           static_solution.load_vector_electric, **kwargs)
        elif isinstance(u0, CurrentFlowSolutionAxiStatic):
            static_solution = u0
        else:
            static_solution = CurrentFlowSolutionAxiStatic.solution_from_problem(self.description, self, u0)
        self.boundary_conditions.update('solution', static_solution)
        # endregion

        # region Newton algorithm
        while True:
            iter_no += 1
            loss_power_prev = static_solution.joule_loss_power

            # region Build and solve system
            # ToDo: Split sigmad in sigma and 2 dsigma/dE^2 --> avoid loop over elements with linear materials
            rhs = coo_matrix(static_solution.load_vector_electric
                             - static_solution.divgrad_matrix_sigma @ static_solution.potential[:, None]
                             + static_solution.divgrad_matrix_sigmad @ static_solution.potential[:, None])
            lhs = static_solution.divgrad_matrix_sigmad
            potential_new = self._solve_system(lhs, rhs, **kwargs)
            # endregion

            # region Insert relaxation and write new potential vector to solution
            if solver_info.relaxation_newton != 1:
                potential_new = solver_info.relaxation_newton * potential_new + \
                    (1 - solver_info.relaxation_newton) * static_solution.potential

            # Update solution and nonlinear quantities
            static_solution.potential = potential_new
            self.boundary_conditions.update('solution', static_solution)
            # endregion

            # region Terminate algorithm
            # Compute relative error
            relative_error = (abs((loss_power_prev - static_solution.joule_loss_power) /
                                  static_solution.joule_loss_power))

            logger.info('Iteration: %d, rel. error: %e, power loss: %f', iter_no, relative_error,
                        static_solution.joule_loss_power)

            # Check for convergence
            if relative_error <= solver_info.tolerance_newton:
                logger.info('Convergence after %d iterations. Relative error: %e', iter_no, relative_error)
                return static_solution

            # Check if maximum number of iterations is exceeded.
            if iter_no == solver_info.max_iter_newton:
                logger.critical('Maximum number of iterations %d exceeded. Relative error: %f', iter_no, relative_error)
                return static_solution
            # endregion
        # endregion
        # endregion


class CurrentFlowProblemAxiHarmonic(HarmonicProblem):
    r"""A harmonic electroquasistatic problem in axisymmetric coordinates:

    The harmonic electroquasistatic problem models capacitive-resistive effects. The corresponding differential \
    equation reads

    .. math::
        -\mathrm{div}(j 2\pi f\varepsilon\, \mathrm{grad} (\phi)) -\mathrm{div}(\sigma\,\mathrm{grad}
         (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`~pyrit.material.Conductivity`),
    :math:`\varepsilon` is the electric permittivity (see :py:class:`~pyrit.material.Permittivity`), :math:`\phi`
    denotes the electric scalar potential, and :math:`f` is the frequency. The problem can, by definition, only handle
    linear materials. A possible application is, for example, the steady state simulation of the insulation system in
    an electrical machine.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.CurrentFlowSolutionAxiHarmonic`.
    """

    problem_identifier: str = 'Harmonic electroquasistatic problem in axisymmetric coordinates'

    def solve(self, *args, **kwargs) -> CurrentFlowSolutionAxiHarmonic:
        raise NotImplementedError()


class CurrentFlowProblemAxiTransient(TransientProblem):
    r"""A transient electroquasistatic problem in axisymmetric coordinates:

    The transient electroquasistatic problem models capacitive-resistive effects. The corresponding differential \
    equation reads

    .. math::
        -\mathrm{div}(\partial_t(\varepsilon\, \mathrm{grad} (\phi))) -\mathrm{div}(\sigma\,\mathrm{grad}
         (\phi)) = 0,

    where :math:`\sigma` is the electric conductivity (see :py:class:`Conductivity`), :math:`\varepsilon` is the \
    electric permittivity (see :py:class:`Permittivity`) and :math:`\phi` denotes the electric scalar potential. \
    A possible application is, for example, the simulation of transient effects in HVDC cable joints or AC surge \
    arresters.

    In case of nonlinear problems, the solve routine performs a Newton algorithm. The algorithm requires the definition
    of a :py:class:`~pyrit.material.DifferentialConductiviy` and :py:class:`~pyrit.material.DifferentialPermittivity`
    for materials with a nonlinear conductivity or nonlinear permittivity, respectively.

    The corresponding solution class is :py:class:`~pyrit.problem.solutions.CurrentFlowSolutionAxiTransient`.
    """

    problem_identifier: str = 'Transient electroquasistatic problem in axisymmetric coordinates'
    available_monitors: dict = {'joule_loss_power': '_joule_loss_power'}

    def __init__(self, description: str, axi_mesh: mesh.AxiMesh,
                 tri_axi_nodal_shape_function: shapefunction.TriAxisymmetricNodalShapeFunction,
                 regions: region.Regions, materials: material.Materials,
                 boundary_conditions: bdrycond.BdryCond, excitations: excitation.Excitations, time_steps: np.ndarray):
        """Transient EQS problem, axisymmetric.

        Parameters
        ----------
        description : str
            A description of the problem.
        axi_mesh : AxiMesh
            A mesh object. See :py:class:`pyrit.mesh.AxiMesh`.
        tri_axi_nodal_shape_function : TriAxisymmetricNodalShapeFunction
            A shape function object. See :py:class:`pyrit.shapefunction.TriAxisymmetricNodalShapeFunction`
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        boundary_conditions : BdryCond
            A boundary conditions object. See :py:mod:`pyrit.bdrycond`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        time_steps : np.ndarray
            An array with the time steps for the problem.
        """
        super().__init__(description, None, None, regions, materials, boundary_conditions, excitations, time_steps)

        self.mesh: mesh.AxiMesh = axi_mesh
        self.shape_function: shapefunction.TriAxisymmetricNodalShapeFunction = tri_axi_nodal_shape_function

        self.divgrad_eps_linear = self.materials.is_linear(material.Permittivity)
        self.divgrad_sigma_linear = self.materials.is_linear(material.Conductivity)
        self.divgrad_epsd_linear = self.materials.is_linear(material.DifferentialPermittivity)
        self.divgrad_sigmad_linear = self.materials.is_linear(material.DifferentialConductivity)
        self.load_linear = self.excitations.is_linear
        self.boundary_condition_linear = self.boundary_conditions.is_linear
        self.all_linear = all([self.divgrad_eps_linear, self.divgrad_epsd_linear, self.divgrad_sigma_linear,
                               self.divgrad_sigmad_linear, self.load_linear, self.boundary_condition_linear])
        self._is_temperature_dependent: bool = None

        self.consistency_check()

    @property
    def is_temperature_dependent(self):
        """True if one of the materials has a conductivity with the key word argument temperature."""
        if self._is_temperature_dependent is None:
            self._is_temperature_dependent = False
            for mat in self.materials:
                if 'temperature' in mat.get_property(material.Conductivity).keyword_args:
                    self._is_temperature_dependent = True
            for exci in self.excitations:
                if 'temperature' in exci.keyword_args:
                    self._is_temperature_dependent = True
            for bc in self.boundary_conditions:
                if 'temperature' in bc.keyword_args:
                    self._is_temperature_dependent = True
        return self._is_temperature_dependent

    def consistency_check(self):
        # make sure all materials with a nonlinear conductivity also have a differential conductivity
        if not self.divgrad_sigma_linear:
            for mat_id in self.materials.get_ids():
                mat = self.materials.get_material(mat_id)
                cond = mat.get_property(material.Conductivity)
                if cond is not None:
                    if not cond.is_linear:
                        if mat.get_property(material.DifferentialConductivity) is None:
                            raise ValueError("All materials with a nonlinear Conductivity must have a "
                                             "DifferentialConductivity. The DifferentialConductivity of "
                                             + mat.name + " is missing.")

        if not self.divgrad_eps_linear:
            for mat_id in self.materials.get_ids():
                mat = self.materials.get_material(mat_id)
                eps = mat.get_property(material.Permittivity)
                if eps is not None:
                    if not eps.is_linear:
                        if mat.get_property(material.DifferentialPermittivity) is None:
                            raise ValueError("All materials with a nonlinear Permittivity must have a "
                                             "DifferentialPermittivity. The DifferentialPermittivity of "
                                             + mat.name + " is missing.")

    @staticmethod
    def _joule_loss_power(solution: 'CurrentFlowSolutionAxiStatic'):
        return solution.joule_loss_power

    Monitor = Union[str, Callable[['CurrentFlowSolutionAxiStatic'], Any]]
    Solution_Monitor = Union[int, Iterable[int]]

    def _solve_system(self, matrix, rhs, **kwargs) -> np.ndarray:
        matrix_shrink, rhs_shrink, _, _, support_data = self.shape_function.shrink(matrix, rhs, self, 1)
        potential_shrink, _ = type(self).solve_linear_system(matrix_shrink.tocsr(), rhs_shrink.todense(),
                                                             **kwargs)

        prev_sol = self.shape_function.inflate(potential_shrink, self, support_data)
        return prev_sol

    def solve(self, start_value: np.ndarray, solution_monitor: Union[int, Iterable[int]] = 1,
              monitors: Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]] = None,
              callback: Callable[['CurrentFlowSolutionAxiStatic'], NoReturn] = None,
              solver_info: CurrentFlowSolverInfoAxiTransient = None, **kwargs) -> \
            CurrentFlowSolutionAxiTransient:

        # region Initialization and Preprocessing
        if not self.all_linear:
            if solver_info is None:
                solver_info = CurrentFlowSolverInfoAxiTransient()

        if not isinstance(start_value, np.ndarray):
            raise ValueError(f"The argument start_value is not of the right type. "
                             f"It is '{type(start_value)}' but should be a ndarray")

        if isinstance(solution_monitor, int):
            solution_monitor = np.arange(start=0, stop=len(self.time_steps), step=solution_monitor)
            if solution_monitor[-1] != len(self.time_steps) - 1:
                solution_monitor = np.concatenate([solution_monitor, np.array([len(self.time_steps) - 1])])

        time_steps = self.time_steps[solution_monitor]
        delta_ts = np.diff(self.time_steps)

        solution_monitor = set(solution_monitor)

        # Preprocess the monitors so that they have a fixed form
        monitors = self._monitors_preprocessing(monitors)
        monitors_results = {k: [self.time_steps[v[0]], []] for k, v in monitors.items()}

        # Add the start value to the solutions
        if 0 in solution_monitor:
            solutions = [start_value, ]
        else:
            solutions = []

        has_callback = callback is not None
        # endregion

        # Loop over all time steps
        for k in range(1, len(self.time_steps)):
            logger.info('starting with time step %d', k)
            delta_t = delta_ts[k - 1]

            # region Initialize static solution
            if k == 1:
                if isinstance(start_value, CurrentFlowSolutionAxiIntermediate):
                    static_solution = start_value
                elif isinstance(start_value, CurrentFlowSolutionAxiStatic):
                    static_solution = CurrentFlowSolutionAxiIntermediate.extend_solution(start_value)
                else:
                    static_solution = CurrentFlowSolutionAxiIntermediate(f'temporary static solution generated by '
                                                                           f'solve in ElectroquasistaticProblemAxi '
                                                                           f'at time step {k}',
                                                                           start_value, self.mesh, self.shape_function,
                                                                           self.regions,
                                                                           self.materials, self.excitations)
                static_solution.time = self.time_steps[k - 1]
                self.boundary_conditions.update('solution', static_solution)
            # endregion

            # Compute part of the rhs that depends only on the previous time step
            rhs_helper = 1 / delta_t * static_solution.divgrad_matrix_eps @ static_solution.potential[:, None]

            # Update materials and matrices to current time step
            static_solution.time = self.time_steps[k]
            self.boundary_conditions.update('time', self.time_steps[k])

            # Save monitor results of first time step
            if k == 1:
                for key, monitor in monitors.items():
                    if k in monitor[0]:
                        monitors_results[key][1].append(monitor[1](static_solution))

            if self.all_linear:
                # region Solve system in case of a linear problem
                # Finish building rhs and system matrix
                rhs = rhs_helper + static_solution.load_vector_electric
                matrix = 1 / delta_t * static_solution.divgrad_matrix_eps + static_solution.divgrad_matrix_sigma
                # Solve system
                static_solution.potential = self._solve_system(matrix, coo_matrix(rhs), **kwargs)
                # endregion

            else:
                # region Perform Newton algorithm in case of a nonlinear problem
                iter_no = 0
                while True:
                    # Update iteration variables
                    iter_no += 1
                    loss_power_prev = static_solution.joule_loss_power
                    # region Build and solve system
                    # ToDo: Split differential permittivity/ conductivity into two different matrices
                    #  (sigma + tensor term) --> loop over fewer elements for Ksigmad
                    rhs = rhs_helper + \
                          static_solution.load_vector_electric - \
                          static_solution.divgrad_matrix_sigma @ static_solution.potential[:, None] + \
                          static_solution.divgrad_matrix_sigmad @ static_solution.potential[:, None] - \
                          1 / delta_t * static_solution.divgrad_matrix_eps @ static_solution.potential[:, None] + \
                          1 / delta_t * static_solution.divgrad_matrix_epsd @ static_solution.potential[:, None]
                    matrix = static_solution.divgrad_matrix_sigmad + 1 / delta_t * static_solution.divgrad_matrix_epsd
                    potential_new = self._solve_system(matrix, coo_matrix(rhs), **kwargs)
                    # endregion

                    # region Insert relaxation and write new potential vector to solution
                    if solver_info.relaxation_newton != 1:
                        potential_new = solver_info.relaxation_newton * potential_new + \
                            (1 - solver_info.relaxation_newton) * static_solution.potential

                    static_solution.potential = potential_new
                    self.boundary_conditions.update('solution', static_solution)
                    # endregion

                    # region Terminate Newton algorithm
                    relative_error = (abs((loss_power_prev - static_solution.joule_loss_power) /
                                          static_solution.joule_loss_power))

                    # Report progress
                    # logger.info('Newton iteration: %d, rel. error: %e, loss power: %f', iter_no,
                    #            relative_error, static_solution.loss_power)

                    # Check for convergence
                    if relative_error <= solver_info.tolerance_newton:
                        logger.info('Convergence after %d Newton iterations. Relative error: %f',
                                    iter_no, relative_error)
                        break

                    # Check if maximum number of iterations is exceeded.
                    if iter_no == solver_info.max_iter_newton:
                        logger.critical('Maximum number of iterations %d exceeded. Relative error: %f', iter_no,
                                        relative_error)
                        break
                    # endregion

                # endregion

            # region Post-process results of current time step
            if has_callback:
                callback(static_solution)

            for key, monitor in monitors.items():
                if k in monitor[0]:
                    monitors_results[key][1].append(monitor[1](static_solution))

            if k in solution_monitor:
                solutions.append(static_solution.potential)
            # endregion

        # region Write output
        solutions = np.stack(solutions)

        for key, monitor in monitors_results.items():
            monitors_results[key][1] = np.stack(monitors_results[key][1])
        # endregion
        return CurrentFlowSolutionAxiTransient(f'solution of {self.description}', solutions, self.mesh,
                                               self.shape_function, self.regions, self.materials,
                                               self.excitations,
                                               time_steps, monitors_results)
