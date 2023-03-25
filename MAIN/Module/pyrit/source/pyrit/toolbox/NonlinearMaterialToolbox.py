# coding=utf-8
"""NonlinearMaterialToolbox

.. sectionauthor:: Ruppert
"""
from typing import Callable, Union
import numpy as np
from pyrit.solution import StaticSolution


class NonlinearProperty:
    """Class implementing a nonlinear material property as a callable object (see :py:class:`pyrit.ValueHandler.py`)."""

    def __init__(self, f: Callable[[StaticSolution, np.ndarray], Union[float, np.ndarray]],
                 solution: StaticSolution = None, element: np.ndarray = 0):
        r"""Constructor of NonlinearProperty

        Parameters
        ----------
        f : Callable[ [StaticSolution, np.ndarray]
            f(solution, element): Function that defines the nonlinear material property.
        solution: StaticSolution
            StaticSolution defining the current working point of the material.
        element: np.ndarray
            (E,)-array of element indices where the value of the property should be evaluated.
        """
        self.solution = solution
        self.element = element
        self.f = f

    def __call__(self, solution: StaticSolution = None, element: Union[np.ndarray, int] = None):
        r"""Callable of NonlinearProperty.

        Parameters
        ----------
        solution: StaticSolution
            StaticSolution defining the current working point of the material.
        element: np.ndarray
            (E,)-array of element indices where the value of the property should be evaluated.

        Returns
        -------
        out: np.ndarray
            (E,)-array of the property evaluated for the given solution and at the given elements.
        """
        if solution is None:
            solution = self.solution
        if element is None:
            element = self.element
        if isinstance(element, float):
            element = int(element)
        return self.f(solution, element)


def get_differential_tensor(f: np.ndarray, dfdeabs2: np.ndarray, e_field: np.ndarray):
    r"""Computes the differential material tensor for a field-dependent material property, i.e.

    .. math::
        \mathrm{differential\_tensor} = f \cdot \mathrm{unit\_tensor}
        + 2 \cdot \frac{\mathrm{d} f}{\mathrm{d} E^2}(\boldsymbol{E}) \cdot
        \boldsymbol{E}\boldsymbol{E}^{\mathrm{T}}.

    .. table:: Symbols

        ======  =======
        Symbol  Meaning
        ======  =======
        f       function describing the nonlinear material property
        E       electric field strength
        e       Number of elements for which the tensor is computed.
        ======  =======

    Parameters
    ----------
    f : np.ndarray
        (e,)-array of the element-wise material property.
    dfdeabs2: np.ndarray
        (e,)-array of the element-wise derivative of f to (E**2)
    e_field: np.ndarray
        (e,2)-array of the electric field.

    Returns
    -------
    differential_tensor : np.ndarray
        (e,2,2)-tensor of the differential material tensor.
    """
    if isinstance(f, np.ndarray):
        num_elem = len(f)
    else:
        num_elem = 1
    differential_tensor = np.zeros((num_elem, 2, 2))
    differential_tensor[:, 0, 0] = f + 2 * e_field[:, 0] ** 2 * dfdeabs2
    differential_tensor[:, 0, 1] = 2 * e_field[:, 0] * e_field[:, 1] * dfdeabs2
    differential_tensor[:, 1, 0] = differential_tensor[:, 0, 1]
    differential_tensor[:, 1, 1] = f + 2 * e_field[:, 1] ** 2 * dfdeabs2

    return differential_tensor


class NonlinearDifferentialProperty:
    """Class implementing differential material tensors.

    Class that implements the differential material tensor of a field-dependent material as a
    callable object. For more information see :py:class:`pyrit.ValueHandler.py`.
    """

    def __init__(self, f: Callable[[StaticSolution, np.ndarray], Union[float, np.ndarray]],
                 dfdE2: Callable[[StaticSolution, np.ndarray], Union[float, np.ndarray]],
                 solution: StaticSolution = None, element: np.ndarray = 0):
        r"""Constructor of NonlinearDifferentialProperty

        .. table:: Symbols

            ======  =======
            Symbol  Meaning
            ======  =======
            E       electric field strength
            e       Number of elements for which the tensor is computed.
            ======  =======

        Parameters
        ----------
        f : Callable[[StaticSolution, np.ndarray], Union[float, np.ndarray]]
            f(solution, element): Function that defines the nonlinear material property.
        dfdE2 : Callable[[StaticSolution, np.ndarray], Union[float, np.ndarray]]
            dfdE2(solution, element): Function that defines the derivative of the nonlinear material property w.r.t.
            (E ** 2)
        solution: StaticSolution
            StaticSolution defining the current working point of the material.
        element: np.ndarray
            (e,)-array of element indices where the value of the property should be evaluated.
        """
        self.solution = solution
        self.element = element
        self.f = f
        self.dfdE2 = dfdE2

    def __call__(self, element: Union[np.ndarray, int] = None, solution: StaticSolution = None):
        r"""Callable of NonlinearDifferentialProperty

        Parameters
        ----------
        solution: StaticSolution
            StaticSolution defining the current working point of the material.
        element: np.ndarray
            (E,)-array of element indices where the value of the property should be evaluated.

        Returns
        -------
        out: np.ndarray
        (E,2,2)-array of the differential property tensor evaluated for the given solution and at the given elements.
        """
        if solution is None:
            solution = self.solution
        if element is None:
            element = self.element
        if not isinstance(element, np.ndarray):
            element = [int(element)]
        return get_differential_tensor(self.f(solution, element), self.dfdE2(solution, element),
                                       solution.e_field[element, :])
