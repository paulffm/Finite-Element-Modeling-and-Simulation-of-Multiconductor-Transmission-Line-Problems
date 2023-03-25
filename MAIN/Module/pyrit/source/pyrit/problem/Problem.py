# coding=utf-8
"""
Implementation of the Problem class

.. sectionauthor:: Bundschuh
"""

from pathlib import Path
from typing import Union, TYPE_CHECKING, List, Tuple, Dict, Any, Callable, NoReturn, Iterable
from abc import ABC, abstractmethod
import pickle
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as splinalg

from pyrit import get_logger

flag_pypardiso = True
try:
    import pypardiso
except ImportError as e:
    flag_pypardiso = False

    class Pypardiso:
        """Alternative class if pypardiso can't be imported. Only purpose of this class is calling an ImportError."""

        def __init__(self, import_error: ImportError) -> None:
            self.e = import_error

        def __call__(self, *args, **kwargs) -> NoReturn:
            raise self.e

        def spsolve(self, *args, **kwargs) -> NoReturn:
            """Raise Import error if spsolve of pypardiso is called."""
            raise self.e

    pypardiso = Pypardiso(e)

if TYPE_CHECKING:
    from . import StaticSolution, HarmonicSolution, TransientSolution

    from pyrit.mesh import Mesh, TriMesh, TetMesh
    from pyrit.region import Regi
    from pyrit.excitation import Excitations
    from pyrit.material import Materials
    from pyrit.bdrycond import BdryCond
    from pyrit.region import Regions
    from pyrit.shapefunction import ShapeFunction

logger = get_logger(__name__)


class SolverUnknownError(Exception):
    """Custom Error when a solver in the function `solve_linear_system` is not known."""


def solve_linear_system(matrix: sparse.spmatrix, rhs: Union[sparse.spmatrix, np.ndarray],
                        solver: Union[str, Callable[[sparse.csr_matrix, sparse.csr_matrix],
                                                    Tuple[np.ndarray, Dict[str, Any]]]] = None,
                        **kwargs) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Solve a system of linear equations.

    General method for solving a system of linear equations. The matrix and the right-hand-side (rhs) must be given.
    With solver, it can be determined how to solve the system. There are some predefined solvers. It can also be passed
    a function, if you want to provide your own solver. The signature has to be as defined.

    Parameters
    ----------
    matrix : sparse.spmatrix
        The system matrix.
    rhs : Union[sparse.spmatrix, np.ndarray]
        The right-hand-side vector
    solver : Union[str, Callable[[sparse.csr_matrix, sparse.csr_matrix], Tuple[np.ndarray, Dict[str, Any]]]], optional
        The used solver. In the case of a string, use one of the predefined solvers. In the case of a function, make
        sure that the signature is as expected. The function gets the system matrix as first argument and the
        right-hand-side as second argument. It must return the solution as first value and a dictionary with information
        as second value.  Furthermore, the kwargs are passed to the function.
    kwargs :
        kwargs passed to the solver.

    Returns
    -------
    solution : np.ndarray
        The solution of the system of linear equations.
    info : Dict[str, Any]
        A dictionary with information from the solver.
    """
    if not solver:
        solver = kwargs.get('solver', 'spsolve') if flag_pypardiso else kwargs.get('solver', 'spsolve')

    matrix = matrix.tocsr()
    if not isinstance(rhs, (sparse.spmatrix, np.ndarray)):
        raise ValueError("The type of the right hand side is not known.")
    if isinstance(rhs, sparse.spmatrix):
        rhs = rhs.tocsr()

    # Iterative solvers from scipy.sparse.linalg. All have the same basic signature
    iteratiev_solvers = {'bicg': splinalg.bicg, 'bicgstab': splinalg.bicgstab, 'cg': splinalg.cg, 'cgs': splinalg.cgs,
                         'gmres': splinalg.gmres, 'lgmres': splinalg.lgmres, 'minres': splinalg.minres,
                         'qmr': splinalg.qmr, 'gcrotmk': splinalg.gcrotmk, 'tfqmr': splinalg.tfqmr}
    if callable(solver):
        if isinstance(rhs, np.ndarray):
            rhs = sparse.csr_matrix(rhs)
        result, info = solver(matrix, rhs, **kwargs)
    else:
        if solver == 'spsolve':
            result = splinalg.spsolve(matrix, rhs, **kwargs)
            info = {}
        elif solver in iteratiev_solvers:
            if isinstance(rhs, sparse.spmatrix):
                rhs = rhs.toarray()
            result, conv_info = iteratiev_solvers[solver](matrix, rhs, **kwargs)
            info = {'convergence info': conv_info}
        elif solver == 'pardiso':
            if np.issubdtype(matrix.dtype, np.complexfloating) or np.issubdtype(rhs.dtype, np.complexfloating):

                tmp_matrix = sparse.vstack([sparse.hstack([np.real(matrix), -1 * np.imag(matrix)]),
                                            sparse.hstack([np.imag(matrix), np.real(matrix)])])
                tmp_rhs = sparse.vstack([np.real(rhs), np.imag(rhs)])
                tmp_result = pypardiso.spsolve(tmp_matrix, tmp_rhs.toarray(), **kwargs)
                tmp_result = np.reshape(tmp_result, (2, -1))
                result = tmp_result[0] + 1j * tmp_result[1]
            else:
                result = pypardiso.spsolve(matrix, rhs, **kwargs)
            info = {}
        else:
            raise SolverUnknownError(f"Solver '{solver}' not known.")

    return result, info


class Problem:
    """General representation of a problem."""

    #: Identifier for the __repr__ method.
    problem_identifier: str = 'General problem'

    def __init__(self, description: str, mesh: 'Mesh', shape_function: 'ShapeFunction', regions: 'Regions',
                 materials: 'Materials', boundary_conditions: 'BdryCond', excitations: 'Excitations'):
        """A general, non-specified problem.

        Parameters
        ----------
        description : str
            A description of the problem
        mesh : Mesh
            A mesh object. See :py:mod:`pyrit.mesh`.
        shape_function : ShapeFunction
            A shape function object. See :py:mod:`pyrit.shapefunction`.
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        boundary_conditions : BdryCond
            A boundary conditions object. See :py:mod:`pyrit.bdrycond`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        """
        self.description = description
        self.mesh = mesh
        self.shape_function = shape_function
        self.regions = regions
        self.materials = materials
        self.boundary_conditions = boundary_conditions
        self.excitations = excitations

    solve_linear_system = solve_linear_system

    def __repr__(self):
        return f"{self.problem_identifier}: {self.description}"

    def consistency_check(self):
        """Check the consistency of the problem.

        It is checked if for every region there is a related material, boundary condition or excitation in the
        corresponding data structure.
        """
        logger.info("Starting general consistency check on problem '%s'.", self.description)
        for region in self.regions:
            region: 'Regi'
            if region.mat:
                try:
                    self.materials.get_material(region.mat)
                except KeyError:
                    logger.warning("The material with ID %d does not exist in materials.", region.mat)
            if region.bc:
                try:
                    self.boundary_conditions.get_bc(region.bc)
                except KeyError:
                    logger.warning("The boundary condition with ID %d does not exist in boundary_conditions.",
                                   region.bc)
            if region.exci:
                try:
                    self.excitations.get_exci(region.exci)
                except KeyError:
                    logger.warning("The excitation with ID %d does not exist in excitations.", region.exci)
        logger.info("Done with general consistency check on problem '%s'.", self.description)

    def save(self, path: Union[str, 'Path'], ignore_attributes: List[str] = None) -> 'Path':
        """Saves the problem instance.

        Parameters
        ----------
        path : Union[str, Path]
            The path. If a string is given, it is converted to a Path object. The file ending can but must not be given.
            In any case it will be a *.pkl* file. See `pickle doc <https://docs.python.org/3/library/pickle.html>`_ for
            more information.
            If the path contains a non-existing folder, this folder will be created.
        ignore_attributes: List[str], optional
            A list of attributes that are not included in the saved file. By default, all attributes are included.

        Returns
        -------
        path : Path
            The path that was used for opening the file.
        """
        path = Path(path).with_suffix('.pkl')

        if not path.parent.exists():
            logger.info("Creating path: %s.", path.parent)
            path.parent.mkdir()

        tmp_dict = {}
        if ignore_attributes:
            for attribute in ignore_attributes:
                try:
                    tmp_dict[attribute] = self.__getattribute__(attribute)
                    delattr(self, attribute)
                except AttributeError:
                    logger.warning("The object does not contain the attribute '%s'. It is ignored.", attribute)

        logger.info("Saving problem to %s.", path)
        with open(path, 'wb') as file:
            pickle.dump(self, file, pickle.HIGHEST_PROTOCOL)
        logger.info("Done saving problem.")

        if ignore_attributes:
            for attribute in ignore_attributes:
                try:
                    self.__setattr__(attribute, tmp_dict[attribute])
                except KeyError:
                    logger.info("The key '%s' does not exist and the object attribute cannot be set.", attribute)

        return path

    @staticmethod
    def load(path: Union[str, Path]) -> 'Problem':
        """Load a problem.

        Parameters
        ----------
        path : Union[str, Path]
            The path. If a string is given, it is converted to a Path object. The file ending can but must not be given.
            In any case it will be a *.pkl* file. See `pickle doc <https://docs.python.org/3/library/pickle.html>`_ for
            more information.

        Returns
        -------
        problem : Problem
            A problem instance.
        """
        path = Path(path).with_suffix('.pkl')

        logger.info("Start reading problem from %s.", path)
        with open(path, 'rb') as file:
            problem = pickle.load(file)
        logger.info("Done reading problem.")
        return problem

    def _check_all_regions_in_trimesh(self):
        """Check if for a TriMesh, every region is included in the mesh."""
        self.mesh: 'TriMesh'
        for regi in self.regions:
            if regi.dim == 0:
                if regi.ID not in self.mesh.node2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)
            elif regi.dim == 1:
                if regi.ID not in self.mesh.edge2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)
            elif regi.dim == 2:
                if regi.ID not in self.mesh.elem2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)

    def _check_all_regions_in_tetmesh(self):
        """Check if for a TetMesh, every region is included in the mesh."""
        self.mesh: 'TetMesh'
        for regi in self.regions:
            if regi.dim == 0:
                if regi.ID not in self.mesh.node2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)
            elif regi.dim == 1:
                if regi.ID not in self.mesh.edge2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)
            elif regi.dim == 2:
                if regi.ID not in self.mesh.face2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)
            elif regi.dim == 3:
                if regi.ID not in self.mesh.elem2regi:
                    logger.warning(
                        "The region with ID %d and dimension %d is not included in the mesh.", regi.ID, regi.dim)


class StaticProblem(Problem, ABC):
    """A general, static problem."""

    problem_identifier: str = 'General static problem'

    def __init__(self, description: str, mesh: 'Mesh', shape_function: 'ShapeFunction', regions: 'Regions',
                 materials: 'Materials', boundary_conditions: 'BdryCond', excitations: 'Excitations'):
        """A general, static problem.

        Parameters
        ----------
        description : str
            A description of the problem
        mesh : Mesh
            A mesh object. See :py:mod:`pyrit.mesh`.
        shape_function : ShapeFunction
            A shape function object. See :py:mod:`pyrit.shapefunction`.
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        boundary_conditions : BdryCond
            A boundary conditions object. See :py:mod:`pyrit.bdrycond`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        """
        super().__init__(description, mesh, shape_function, regions, materials, boundary_conditions, excitations)

    @abstractmethod
    def solve(self, *args, **kwargs) -> 'StaticSolution':
        """Solve the problem"""


class HarmonicProblem(Problem, ABC):
    """A general, harmonic problem."""

    problem_identifier: str = 'General harmonic problem'

    def __init__(self, description: str, mesh: 'Mesh', shape_function: 'ShapeFunction', regions: 'Regions',
                 materials: 'Materials', boundary_conditions: 'BdryCond', excitations: 'Excitations',
                 frequency: float):
        """A general, harmonic problem.

        Parameters
        ----------
        description : str
            A description of the problem
        mesh : Mesh
            A mesh object. See :py:mod:`pyrit.mesh`.
        shape_function : ShapeFunction
            A shape function object. See :py:mod:`pyrit.shapefunction`.
        regions : Regions
            A regions object. See :py:mod:`pyrit.regions`.
        materials : Materials
            A materials object. See :py:mod:`pyrit.materials`.
        boundary_conditions : BdryCond
            A boundary conditions object. See :py:mod:`pyrit.bdrycond`.
        excitations : Excitations
            An excitations object. See :py:mod:`pyrit.excitation`.
        frequency : float
            The frequency of the problem
        """
        super().__init__(description, mesh, shape_function, regions, materials, boundary_conditions, excitations)
        self._frequency = None
        self.frequency = frequency

    @property
    def frequency(self) -> float:
        """The frequency."""
        return self._frequency

    @frequency.setter
    def frequency(self, frequency):
        """The frequency.

        Parameters
        ----------
        frequency : float
            The non-negative frequency.
        """
        if frequency < 0:
            raise ValueError("The frequency must be greater or equal to 0")
        self._frequency = frequency

    @property
    def angular_frequency(self) -> float:
        """The angular frequency."""
        return 2 * np.pi * self.frequency

    omega = angular_frequency

    @abstractmethod
    def solve(self, *args, **kwargs) -> 'HarmonicSolution':
        """Solve the problem"""


Monitor = Union[str, Callable[['StaticSolution'], Any]]
Solution_Monitor = Union[int, Iterable[int]]  # pylint: disable=invalid-name


class TransientProblem(Problem, ABC):
    """A general, transient problem."""

    problem_identifier: str = 'General transient problem'
    available_monitors: dict = {}

    def __init__(self, description: str, mesh: 'Mesh', shape_function: 'ShapeFunction', regions: 'Regions',
                 materials: 'Materials', boundary_conditions: 'BdryCond', excitations: 'Excitations',
                 time_steps: np.ndarray):
        """A general, harmonic problem.

        Parameters
        ----------
        description : str
            A description of the problem
        mesh : Mesh
            A mesh object. See :py:mod:`pyrit.mesh`.
        shape_function : ShapeFunction
            A shape function object. See :py:mod:`pyrit.shapefunction`.
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
        super().__init__(description, mesh, shape_function, regions, materials, boundary_conditions, excitations)
        self.time_steps = time_steps

    @property
    def start_time(self):
        """The start time of the problem."""
        return self.time_steps[0]

    @property
    def end_time(self):
        """The end time of the problem."""
        return self.time_steps[-1]

    @property
    def time_span(self):
        """The time span of the problem."""
        return self.end_time - self.start_time

    def _monitors_preprocessing(self, monitors):
        if monitors is None:
            return {}

        internal_monitors = {}
        for key, monitor in monitors.items():
            tmp_monitor = [None, None]
            if isinstance(monitor, (list, tuple)):
                if isinstance(monitor[0], int):
                    tmp = np.arange(start=0, stop=len(self.time_steps), step=monitor[0])
                    # solution_monitor = {n for n in range(0, len(self.time_steps), solution_monitor)}
                    # if len(self.time_steps)-1 not in solution_monitor:
                    #     solution_monitor.add(len(self.time_steps)-1)
                    if tmp[-1] != len(self.time_steps) - 1:
                        tmp = np.concatenate([tmp, np.array([len(self.time_steps) - 1])])
                    # tmp = {n for n in range(0, len(self.time_steps), monitor[0])}
                    # if len(self.time_steps) - 1 not in tmp:
                    #     tmp.add(len(self.time_steps) - 1)
                    tmp_monitor[0] = tmp
                else:
                    tmp_monitor[0] = monitor[0]
                monitor = monitor[1]
            else:
                tmp_monitor[0] = np.arange(len(self.time_steps))

            if isinstance(monitor, str):
                tmp_monitor[1] = self.__getattribute__(self.available_monitors[monitor])
            else:
                tmp_monitor[1] = monitor

            internal_monitors[key] = tmp_monitor

        return internal_monitors

    # todo: Here we need more information for the handling of nonlinear problems: explain the Newton and successive
    #  substitution handling or link tutorial examples.

    @abstractmethod
    def solve(self, start_value: np.ndarray, solution_monitor: Solution_Monitor = 1,
              monitors: Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]] = None,
              callback: Callable[['StaticSolution'], NoReturn] = None, **kwargs) -> 'TransientSolution':
        """Solve the problem.

        The problem is defined in the object. A start value has to be given. Furthermore, you can determine at which
        time steps the solution should be stored in the solution object (`solution_monitor`). With `monitors` you can
        determine what entities over time should be calculated during the simulation. Therefore, the solution in each
        time step is used. This can save memory.

        Parameters
        ----------
        start_value : np.ndarray
            The start value.
        solution_monitor : Union[int, np.ndarray]
            Determines at which time steps the solution should be stored in the solution object. If it is an integer
            :math:`n`, the solution is stored at every :math:`n`-th time step (including the first and the last one). If
            it is an array, store at the indicated time steps.
        monitors : Dict['str', Union[Monitor, Tuple[Solution_Monitor, Monitor]]]
            A number of monitors. Each monitor, i.e. each entry in the dict, has a name (the key). The value is either a
            'Monitor' or a tuple (or a list) with the first entry containing the time steps when the monitor is
            evaluated and the second entry being th 'Monitor'. Here, a 'Monitor' is eiter the name of a predefined
            monitor or a function that returns the information. The predefined monitors are saved in the attribute
            `available_monitors`. See the examples for more information.
        callback : Callable[['Solution'], NoReturn]
            If given, this function is executed after every iteration. This can be used to update information in e.g.
            material data.
        kwargs :
            Additional keyword arguments passed to the py:func:`solve_lgs` function.

        Returns
        -------
        solution : TransientSolution
            The solution of the problem.

        Examples
        --------
        Suppose the transient problem is defined in `problem`. The problem shall be simulated in the interval from 0 to
        1 second:

        >>> problem.time_steps = np.linspace(0,1,100)
        >>> start_value = np.zeros(problem.mesh.num_node)

        However, for further computations the solution vector (which will be computed in every of this 100 time steps)
        is only needed at the last 10 time steps. Then, the function call would be:

        >>> solution = problem.solve(start_value, solution_monitor = np.arange(90,100))

        When one wants to have the energy over time, but still save the solution only at the last 10 time steps, the
        function call would be:

        >>> solution = problem.solve(start_value, solution_monitor = np.arange(90,100), {'energy': 'energy'})

        Here, it is assumed that on the problem class there is a predefined monitor called 'energy' defined. If this is
        not the case, one can simply define this monitor by oneself:

        >>> def energy_monitor(static_solution):
        >>>     return static_solution.energy()
        >>>
        >>> solution = problem.solve(start_value, solution_monitor = np.arange(90,100), {'energy': energy_monitor})

        """
