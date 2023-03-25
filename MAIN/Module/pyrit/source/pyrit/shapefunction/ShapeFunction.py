# -*- coding: utf-8 -*-
"""Contains the abstract class ShapeFunction.

.. sectionauthor:: Christ, Bundschuh
"""

from typing import Any, Tuple, TYPE_CHECKING, Protocol, Dict, Union, Callable, Type
from inspect import signature, isclass
from abc import ABC, abstractmethod
from functools import partial

from numpy import ndarray
import numpy as np
from scipy.sparse import coo_matrix

from pyrit.material import MatProperty, Materials
from pyrit.excitation import Excitations
from pyrit.region import Regions
from pyrit.bdrycond import BdryCond
from pyrit.toolbox.QuadratureToolbox import gauss_triangle, gauss_tetrahedron

if TYPE_CHECKING:
    from pyrit.mesh import Mesh


# pylint: disable=no-else-return


def is_valid_index_array(x: Any) -> bool:
    """Returns true if x is an array with positive integers"""
    if isinstance(x, np.ndarray):
        if np.ndim(x) > 0:
            if issubclass(x.dtype.type, np.integer) and np.min(x) >= 0:
                return True
    return False


def is_single_float(x: Any) -> bool:
    """Returns true if x is an int, float or a one-dimensional numpy array."""
    if isinstance(x, (int, float, complex, np.int64, np.int32)):
        return True
    if isinstance(x, ndarray):
        if len(x.shape) == 1:
            return True
    return False


class ShrinkInflateProblemShapeFunction(Protocol):
    """A protocol for the minimal content of a problem such that the shrink and inflate methods can work"""

    mesh: 'Mesh'
    regions: 'Regions'
    boundary_conditions: 'BdryCond'


class ShapeFunction(ABC):
    """Abstract class for shape function (sft) object.

    See Also
    --------
    pyrit.mesh : Abstract underlying Mesh object on whose entities the SFT
                  parameters are computed.
    """

    @abstractmethod
    def __init__(self, mesh: 'Mesh', dim: int, allocation_size: int) -> None:
        """
        Constructor for the abstract class ShapeFunction.

        Parameters
        ----------
        mesh : Mesh object
            Representing the underlying geometry on which the SFTs are defined.
        dim : int
            Dimensionality of the shape function.
        allocation_size : int
            Size used for vector allocation in matrix creation.

        Returns
        -------
        None
        """
        self.mappingcoeff = None
        self._mesh = mesh
        self._dim = dim
        self.__allocation_size = allocation_size

    @property
    def mesh(self) -> 'Mesh':
        """Returns the mesh."""
        return self._mesh

    @property
    def dim(self) -> int:
        """Return the dimension."""
        return self._dim

    @property
    def _allocation_size(self):
        """Get the size need to allocate for sparse creation."""
        return self.__allocation_size

    def _process_arguments(self, args: Any, length: int = None) -> str:
        """
        Checks the type and the right format of the arguments.

        The expected form of args is for length=3:
            Union[Regions, Materials, MatProperty]]
        and for length=2:
            Union[Regions, Excitations]]

        Parameters
        ----------
        args : Any
            The input argument to check
        length : int, optional
            The expected length of args if args is a tuple

        Returns
        -------
        case : str
            Descriptive string for the input argument, one of {'function', 'array', 'number', 'tuple'}

        Raises
        ------
        ValueError
            When args has not the expected format and type
        """
        if length not in (2, 3, None):
            raise ValueError(f"Variable 'length' has value {length} but has to be 2, 3 or empty.")
        if callable(args):
            if len(signature(args).parameters) != self.mesh.node.shape[1]:
                raise ValueError(f"The function has not the expected number of {str(self.mesh.node.shape[1])} "
                                 f"arguments")
            if is_single_float(args(*self.mesh.node[0])):  # returns a single floating number
                return "function"
            raise ValueError("The function does not return scalar fields (i.e. single, number valued return value).")

        if isinstance(args, ndarray):
            if len(args.shape) != 1:
                if len(args.shape) != 2:
                    raise ValueError("The array has not the right shape")
                if args.shape[1] != 1:
                    raise ValueError("The array has not the right shape")
            if args.shape[0] != self.mesh.num_elem:
                raise ValueError(f"The array has {str(args.shape[0])} entries, but needs {str(self.mesh.num_elem)} "
                                 f"entries")
            return "array"
        if isinstance(args, (float, int, complex)):
            return "number"
        if isinstance(args, tuple):
            if length is None:
                raise ValueError("No length parameter provided to function call. Either the input should not be a "
                                 "tuple or the testing is called wrongly")
            if len(args) != length:
                raise ValueError(f"Wrong number of values. It is {str(len(args))} but should be {str(length)}.")
            if not isinstance(args[0], Regions):
                raise ValueError(f"The first value is of type {str(type(args[0]))} but should be type Regions")
            if length == 3:
                if not isinstance(args[1], Materials):
                    raise ValueError(f"The second value is of type {str(type(args[1]))} but should be type Materials")
                if isclass(args[2]):
                    try:
                        tmp = args[2].__call__(0)
                    except Exception as e:
                        raise ValueError("The third value is not a subclass of MatProperty") from e
                    if not isinstance(tmp, MatProperty):
                        raise ValueError(f"The third value is of type {str(type(tmp))} but should be  a subclass of "
                                         f"MatProperty")
                else:
                    raise ValueError("The third value is not a class")
            elif length == 2:
                if not isinstance(args[1], Excitations):
                    raise ValueError(f"The second value is of type {str(type(args[1]))} but should be type Excitations")
            return "tuple"
        raise ValueError("The argument does not match any expected value")

    def _process_single_argument(self, arg: Any) -> str:
        """Process the arguments if the arguments consists of a single value.

        Parameters
        ----------
        arg : Any
            The single argument.

        Returns
        -------
        type : str
            The type of the argument. Either "function", "array" or "number".

        Raises
        ------
        ValueError
            When arg has not the expected format and type
        """
        if callable(arg):
            if len(signature(arg).parameters) != self.mesh.node.shape[1]:
                raise ValueError(f"The function has not the expected number of {str(self.mesh.node.shape[1])} "
                                 f"arguments")
            if is_single_float(arg(*self.mesh.node[0])):  # returns a single floating number
                return "function"
            raise ValueError("The function does not return scalar fields (i.e. single, number valued return value).")

        if isinstance(arg, ndarray):
            if len(arg.shape) == 1:
                if arg.shape[0] == self.mesh.num_elem:
                    return "array"
                else:
                    raise ValueError(f"The array has {str(arg.shape[0])} entries, but needs {str(self.mesh.num_elem)} "
                                     f"entries")
            if len(arg.shape) == 2:
                if arg.shape[0] == self.mesh.num_elem and arg.shape[1] == 1:
                    return "array"
                if arg.shape[0] == 2 and arg.shape[1] == 2:
                    return "single_tensor"
                else:
                    raise ValueError("The array has not the right shape")
            if len(arg.shape) == 3:
                if arg.shape[0] == self.mesh.num_elem and arg.shape[1] == 2 and arg.shape[2] == 2:
                    return "tensor_per_element"
                else:
                    raise ValueError("The array has not the right shape")

        if isinstance(arg, (float, int)):
            return "number"

        raise ValueError("Argument type not excepted.")

    def _process_material(self, *material: Any) -> Tuple[str, Union[Callable[..., float], ndarray, float,
                                                                    Union[Regions, Materials, Type[MatProperty]]]]:
        """Process the material arguments.

        Parameters
        ----------
        material : Any
            The material arguments.

        Returns
        -------
        type : str
            The type of the arguments. Either "tuple", "function", "array" or "number".
        material : Tuple[str, Union[Callable[..., float], ndarray, float, Union[Regions, Materials, Type[MatProperty]]]]
            The arguments in a standard format. If case is "tuple", the order in `material` is as given (independent of
            the order of the input). If case is not "tuple", `material` is no list but just the value.

        Raises
        ------
        ValueError
            When args has not the expected format and type
        """
        if len(material) not in (1, 3):
            raise ValueError(f"Wrong number of inputs. Got {len(material)} but has to be 1 or 3.")
        material_out = None
        case = ''

        if isinstance(material[0], tuple):  # needed for TetCartesianNodal
            material = material[0]

        if len(material) == 1:
            case = self._process_single_argument(*material)  # pylint: disable=no-value-for-parameter
            material_out = material[0]

        if len(material) == 3:
            case = 'tuple'
            material_out = [None, None, None]
            classes: list = [Regions, Materials, MatProperty]
            for k, mat in enumerate(material):
                if isinstance(mat, Regions):
                    material_out[0] = mat
                elif isinstance(mat, Materials):
                    material_out[1] = mat
                elif isclass(mat):
                    if issubclass(mat, MatProperty):
                        material_out[2] = mat
                else:
                    raise ValueError(f"Unexpected argument: {mat}")

            for k, c in enumerate(classes):
                if material_out[k] is None:
                    raise ValueError(f"The class {c} is missing in the materials.")

        return case, material_out

    def _process_load(self, *load: Any) -> Tuple[str, Union[Callable[..., float], ndarray, float,
                                                            Union[Regions, Excitations]]]:
        """Process the load arguments.

        Parameters
        ----------
        load : Any
            The load arguments.

        Returns
        -------
        type : str
            The type of the arguments. Either "tuple", "function", "array" or "number".
        load : Tuple[str, Union[Callable[..., float], ndarray, float, Union[Regions, Excitations]]]
            The arguments in a standard format. If case is "tuple", the order in `load` is as given (independent of
            the order of the input). If case is not "tuple", `load` is no list but just the value.

        Raises
        ------
        ValueError
            When args has not the expected format and type
        """
        if len(load) not in (1, 2):
            raise ValueError(f"Wrong number of inputs. Got {len(load)} but has to be 1 or 2.")

        load_out = None
        case = ''
        if len(load) == 1:
            case = self._process_single_argument(*load)
            load_out = load[0]

        if len(load) == 2:
            case = 'tuple'
            load_out = [None, None]
            classes: list = [Regions, Excitations]
            for k, l in enumerate(load):
                if isinstance(l, Regions):
                    load_out[0] = l
                elif isinstance(l, Excitations):
                    load_out[1] = l
                else:
                    raise ValueError(f"Unexpected argument: {l}")

            for k, c in enumerate(classes):
                if load_out[k] is None:
                    raise ValueError(f"The class {c} is missing in load.")

        return case, load_out

    def _process_neumann(self, *args: Any, allow_indices_tuple: bool = False) -> Tuple[bool, str]:
        """
        Process the neumann_term arguments.

        Parameters
        ----------
        args: Any
            The neumann_term arguments.
        allow_indices_tuple : bool
            Flag to allow the second argument to be a tuple in case of index value pairs.

        Returns
        -------
        flag_regions : bool
            True, if regions are used.
        flag_value : str
            'none' if regions are used,
            'callable' if indices and a function are passed,
            'array' if indices and a numpy array are passed,
            'value' if indices and a single value are passed.

        Raises
        ------
        ValueError
            If args has not the expected format and type
        Exception
            If len(args) is unequal to two or the first argument is neither a Regions-object nor an array.
        """
        flag_regions = False
        flag_value = 'none'
        if len(args) != 2:
            raise Exception(f"Not right number of arguments given. Needs 2 and has {len(args)}.")
        if isinstance(args[0], Regions) and isinstance(args[1], BdryCond):
            flag_regions = True
        elif isinstance(args[0], ndarray):
            if is_valid_index_array(args[0]) is False:
                raise ValueError("The (first array) provided is not a valid index array.")
            if allow_indices_tuple and isinstance(args[1], tuple):  # used in TriAxisymmetricEdgeShapeFunction
                if len(args[1]) != 1:
                    raise Exception(f"If the second argument is tuple, then its length is supposed to be 1 instead of"
                                    f"{len(args[1])}!")
                arg = args[1][0]
            else:
                arg = args[1]

            if callable(arg):
                num_args = len(signature(arg).parameters)
                if isinstance(arg, partial):
                    # If the provided callable is a functools.partial object, the args and kwargs set in order to
                    # create the partial object are stored in two dicts. The length of the new signature is still the
                    # same as the original function. Hence, num_args that still need to be provided needs to be adapted.
                    num_args -= len(arg.keywords) + len(arg.args)
                if num_args != self.mesh.node.shape[1]:
                    raise ValueError(f"The function has not the expected number of {str(self.mesh.node.shape[1])} "
                                     f"arguments")
                flag_value = 'callable'
            elif isinstance(arg, ndarray):
                if args[0].size != arg.size:
                    raise ValueError(
                        "The value array does not provide exactly one entry per element in the index array.")
                flag_value = 'array'
            elif isinstance(arg, (float, int)):
                flag_value = 'value'
            else:
                raise ValueError(f"Second argument has wrong type (is {type(args[1])})")
        else:
            raise Exception("Arguments have not the required types.")

        return flag_regions, flag_value

    def _default_integrator(self, integration_order: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the weights and local coordinates using the default integrators: gaussian triangle and tetrahedron.

        Parameters
        ----------
        integration_order : int
            Integration Order.

        Returns
        -------
        weights: np.ndarray
            (K,) array. K weights used for numerical integration.
        local_coordinates: np.ndarray
            (K, dim) array. K local coordinates used for numerical integration.
        """
        if self.dim == 2:
            return gauss_triangle(integration_order)
        if self.dim == 3:
            return gauss_tetrahedron(integration_order)
        raise NotImplementedError(f"There is no default integrator implemented for dimension{self.dim}.")

    def _matrix_routine(self, matrix_type: str, *material: Union[Callable[..., float], np.ndarray, float,
                                                                 Union['Regions', 'Materials', 'MatProperty']],
                        integration_order: int = 1, **kwargs) -> coo_matrix:
        """
        Create a matrix according to the matrix assembly routine for a shapefunction-class.

        See Also: doc/src/writing%20a%20shape%20function%20class.html#matrix-assembly-routines.

        Parameters
        ----------
        eval_func : Callable
            Method that handels all the evaluations.
        calc_func : Callable
            Method that handels all the calculations.
        weights : np.ndarray
            Weights for numerical integration.
        local_coordinates : np.ndarray
            Local coordinates for numerical integration.
        material : Union[Callable[..., float], np.ndarray, float, Union['Regions', 'Materials', 'MatProperty']])
            Material.

        Returns
        -------
        coo_matrix
            Matrix in sparse format.
        """
        case, material = self._process_material(*material)

        if not isinstance(integration_order, int):
            raise ValueError(f"integration order is of type {type(integration_order)} but has to be an int.")
        weights, local_coordinates = self._default_integrator(integration_order)  # TO-DO: Consider other integrators

        # allocation for sparse creation:
        i = np.zeros(self._allocation_size, dtype=np.int_)  # row indices
        j = np.zeros(self._allocation_size, dtype=np.int_)  # column indices
        v = np.zeros(self._allocation_size, dtype=np.float_)  # values

        if case == "tuple":
            # using pyrit structures as a material. Thus, multiple regions and different material types are possible.
            regions, materials, mat_property_class = material

            for key in regions.get_keys():  # iterate over regions
                regi = regions.get_regi(key)

                if regi.dim != self.dim:  # only consider regions with the same dimensionality as the shape function
                    continue

                material_id = regi.mat
                if material_id is None:  # only consider regions with a material
                    continue

                material = materials.get_material(material_id)  # get material of current region
                prop = material.get_property(mat_property_class)  # get property of material, e.g. Permeability
                if prop is None:  # Consider only regions that have the desired material property
                    continue

                # get indices of all mesh elements within the region:
                indices = np.where(self.mesh.elem2regi == regi.ID)[0].astype(np.int32, casting="same_kind", copy=False)

                if prop.is_homogeneous:
                    if prop.is_linear:

                        if prop.is_isotrop:  # homogenous, linear, isotropic
                            if prop.is_time_dependent:
                                value = prop.value()
                            else:  # homogenous, linear, isotropic, time independent
                                value = prop.value

                            if isinstance(value, np.ndarray):
                                if len(value.shape) == 1:  # One scalar value per element
                                    self._calc_matrix_scalar_per_elem(matrix_type, indices, i, j, v, value,
                                                                      weights, local_coordinates, **kwargs)
                                elif len(value.shape) == 3:  # One tensor per element
                                    self._calc_matrix_tensor_per_elem(matrix_type, indices, i, j, v, value,
                                                                      weights, local_coordinates, **kwargs)
                                else:
                                    raise NotImplementedError("The value shape is not supported.")
                            else:  # scalar value
                                self._calc_matrix_constant_scalar(matrix_type, indices, i, j, v, value,
                                                                  weights, local_coordinates, **kwargs)

                        else:  # homogenous, linear, anisotropic
                            self._calc_matrix_constant_tensor(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates, evaluator=None, **kwargs)

                    else:  # nonlinear
                        if prop.is_isotrop:  # homogenous, nonlinear, isotropic
                            self._calc_matrix_scalar_per_elem(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_hom_nonlin_iso", **kwargs)
                        else:  # homogenous, nonlinear, anisotropic
                            self._calc_matrix_tensor_per_elem(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_hom_nonlin_aniso", **kwargs)

                else:  # inhomogeneous
                    if prop.is_linear:  # inhomogeneous, linear
                        if prop.is_isotrop:  # inhomogeneous, linear, isotropic
                            self._calc_matrix_function_scalar(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_inhom_lin_iso", **kwargs)
                        else:  # inhomogeneous, linear, anisotropic
                            self._calc_matrix_function_tensor(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_inhom_lin_aniso", **kwargs)
                    else:  # nonlinear
                        if prop.is_isotrop:  # inhomogeneous, nonlinear, isotropic
                            self._calc_matrix_function_scalar(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_inhom_nonlin_iso", **kwargs)
                        else:  # inhomogeneous, nonlinear, anisotropic
                            self._calc_matrix_function_tensor(matrix_type, indices, i, j, v, prop.value,
                                                              weights, local_coordinates,
                                                              evaluator="eval_inhom_nonlin_aniso", **kwargs)

        elif case == "number":  # material is a float -> constant for all mesh elements
            self._calc_matrix_constant_scalar(matrix_type, np.arange(self.mesh.num_elem), i, j, v, float(material),
                                              weights, local_coordinates, **kwargs)
        elif case == "array":  # material is an array -> one value per mesh element
            self._calc_matrix_scalar_per_elem(matrix_type, np.arange(self.mesh.num_elem), i, j, v,
                                              material.astype(float), weights, local_coordinates, **kwargs)
        elif case == "function":  # material is a function
            self._calc_matrix_function_scalar(matrix_type, np.arange(self.mesh.num_elem), i, j, v, material,
                                              weights, local_coordinates, evaluator="eval_inhom_lin_iso", **kwargs)
        elif case == "tensor_per_element":
            self._calc_matrix_tensor_per_elem(matrix_type, np.arange(self.mesh.num_elem), i, j, v, material,
                                              weights, local_coordinates, **kwargs)
        elif case == "single_tensor":
            self._calc_matrix_constant_tensor(matrix_type, np.arange(self.mesh.num_elem), i, j, v, material,
                                              weights, local_coordinates, **kwargs)
        else:
            raise ValueError("Argument type not expected. Method can not handle the case of a " + case)

        return coo_matrix((v, (i, j)), shape=(self.mesh.num_node, self.mesh.num_node))

    @abstractmethod
    def _calc_matrix_constant_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: float,
                                     weights: np.ndarray, local_coordinates: np.ndarray, *args,
                                     evaluator: str = None, **kwargs):
        """
        Calculate a matrix for a material that is represented by a scalar constant.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : float
            Material value.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def _calc_matrix_constant_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray,
                                     weights: np.ndarray, local_coordinates: np.ndarray, *args,
                                     evaluator: str = None, **kwargs) -> None:
        """
        Calculate a matrix for a material that is represented by a constant tensor.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : np.ndarray
            (D,) array. Material value.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def _calc_matrix_scalar_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray,
                                     weights: np.ndarray, local_coordinates: np.ndarray,
                                     *args, evaluator: str = None, **kwargs) -> None:
        """
        Calculate a matrix for a material that is represented by a scalar for each mesh element.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : np.ndarray
            (I,) array. Material value.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def _calc_matrix_tensor_per_elem(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: np.ndarray,
                                     weights: np.ndarray, local_coordinates: np.ndarray, *args,
                                     evaluator: str = None, **kwargs) -> None:
        """
        Calculate a matrix for a material that is represented by a tensor for each mesh element.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : np.ndarray
            (I, D) array. Material value.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def _calc_matrix_function_scalar(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., float],
                                     weights: np.ndarray, local_coordinates: np.ndarray, *args,
                                     evaluator: str = None, **kwargs) -> None:
        """
        Calculate a matrix for a material that is represented by a scalar function.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : Callable
            Material function that returns a scalar value.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def _calc_matrix_function_tensor(self, matrix_type: str, indices: np.ndarray, i: np.ndarray, j: np.ndarray,
                                     v: np.ndarray, value: Callable[..., np.ndarray],
                                     weights: np.ndarray, local_coordinates: np.ndarray, *args,
                                     evaluator: str = None, **kwargs) -> None:
        """
        Calculate a matrix for a material that is represented by a tensor function.

        This method calls the correct calc method for a given matrix type. And passes on the allocated vectors for the
        sparse creation.

        .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        I       Number of indices of mesh elements for which the matrix is computed.
        A       Allocation size.
        D       Dimension of shape function.
        N       Number of evaluation points.
        ======  =======

        Parameters
        ----------
        matrix_type : str
            Label used for identification of matrix type (e.g. "divgard", "curlcurl", "mass").
        indices : np.ndarray
            (I,) array. Indices of mesh elements for which the matrix is computed.
        i : np.ndarray
            (A,) array. Row index vector used for sparse creation.
        j : np.ndarray
            (A,) array. Column index vector used for sparse creation.
        v : np.ndarray
            (A,) array. Value vector used for sparse creation.
        value : Callable
            Material function that returns a tensor.
        weights : np.ndarray
            (N,) array. Weights for each evaluation point (numerical integration).
        local_coordinates : np.ndarray
            (N, D) array. Local coordinates for each evaluation point (numerical integration).
        args : Tuple
            Further arguments passed to the calc method.
        evaluator : str
            Label used to select the correct evaluator method.
        kwargs : Dict
            Further keyword arguments passed to the calc method.

        Returns
        -------
        None
        """

    @abstractmethod
    def shrink(self, matrix: 'coo_matrix', rhs: 'coo_matrix', problem: ShrinkInflateProblemShapeFunction,
               integration_order: int = 1) -> Tuple['coo_matrix', 'coo_matrix', ndarray, int, Dict['str', Any]]:
        r"""Shrink the system of equations of a problem.

        Shrinks the system of equations where `matrix` is the matrix and `rhs` is the right-hand-side. Additionally,
        there is returned side information that can be used for further computations with the system or for inflating.

        Parameters
        ----------
        matrix : coo_matrix
            The matrix of the system of equations.
        rhs : coo_matrix
            The right-hand-side of the system of equations.
        problem : ShrinkInflateProblemShapeFunction
            A problem containing additional data structures. See :py:class:`ShrinkInflateProblem`.
        integration_order : int, optional
            The integration order used by specific shrink methods. Default is 1

        Returns
        -------
        matrix_shrink: coo_matrix
            The shrunk matrix of the system.
        rhs_shrink: coo_matrix
            The shrunk right-hand-side of the system.
        indices_not_dof: ndarray
            An array containing these lines of the original matrix that are no longer in the resulting matrix.
        num_new_dof: int
            The number of newly introduced degrees of freedom
        support_data: Dict[str,Any]
            A dict containing additional data that is needed for inflating the problem.

        Notes
        -----
        Suppose the original matrix :math:`\mathbf{A}` is of dimension :math:`N\times N` and the length of
        `indices_not_dof` is :math:`M`. Then, the dimension of the returned matrix is :math:`N_{out}\times N_{out}` with

        .. math::

            N_{out} = N - M + \text{num\_new\_dof}\,.
        """

    @abstractmethod
    def inflate(self, solution: ndarray, problem: ShrinkInflateProblemShapeFunction,
                support_data: Dict[str, Any] = None) -> ndarray:
        """Inflates a solution.

        After the system of equations of a problem was shrunk and solved, its solution can be inflated with this method.

        Parameters
        ----------
        solution : ndarray
            The solution of the shrunk system of equations.
        problem : ShrinkInflateProblemShapeFunction
            A problem containing additional data structures. See :py:class:`ShrinkInflateProblem`.
        support_data : Dict[str,Any], optional
            A dictionary with data needed for the inflation process. It is provided by the shrink method. If not given,
            these data is computed inside the function. By passing the data provided by the shrink method, this time
            can be saved. Default is None

        Returns
        -------
        solution_inflated: ndarray
            The solution of the original system of equations.
        """
