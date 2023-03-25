# coding=utf-8
"""Class ValueHandler

.. sectionauthor:: Bundschuh
"""

from abc import ABC, abstractmethod
from typing import Any
from inspect import signature

from .Logger import get_logger

logger = get_logger(__name__)


class ValueHandler(ABC):
    """Abstract class that manages a value.

    Values can be numbers, arrays or functions. Depending on the type of values and, if it is a function, on the
    signature, it is interpreted as homogeneous, linear or time dependent. The decision is made as follows:

        - **homogeneous**: It is inhomogeneous only if it is a function with positional arguments (no default)
        - **linear**: If is nonlinear only if it is a function with keyword arguments (has a default) except the
          keyword argument `time`
        - **time dependent**: It is time dependent only if it is a function with the keyword argument `time`.

    About the positional and keyword arguments: It is not distinguished between positional-only or keyword-only
    arguments. In ths scope of this class, an argument is positional if it has no default and a keyword argument if it
    has a default.

    If the value is nonlinear or time dependent, it is strongly recommended to implement the value as a callable object
    of a class. Only then, the method `update` can be used. Furthermore, it is recommended to have an argument named
    'time' (in case of a time dependent value) and an argument named 'solution' (in case of a nonlinear value). So, it
    can be utilized by the solve-routines.

    Examples
    --------
    An example for a value class is

    >>> class Value:
    >>>     def __init__(self, time=0, solution=None):
    >>>         self.time = time
    >>>         self.solution = solution
    >>>
    >>>     def __call__(self, x, y, time=None, solution=None):
    >>>         if time is None:
    >>>             time = self.time
    >>>         if solution is None:
    >>>             solution = self.solution
    >>>         ...

    This value is nonlinear (because of the `solution` parameter), time dependent (because of the `time` parameter) and
    inhomogeneous (because of the positional arguments `x` and `y`).
    """

    __slots__ = ('_value', 'positional_args', 'keyword_args')

    @abstractmethod
    def __init__(self, value):
        self._value = None
        self.positional_args = []
        self.keyword_args = []

        self.value = value

    @staticmethod
    def compute_args(value):
        """Compute the lists of arguments for a given value."""
        positional_args = []
        keyword_args = []

        if callable(value):
            for key, val in signature(value).parameters.items():
                if val.default is val.empty:
                    positional_args.append(key)
                else:
                    keyword_args.append(key)

        if 't' in keyword_args:
            logger.warning("There is a keyword argument 't' in the value. Note that this is not interpreted as time. "
                           "If you wish to have a time dependency, you have to name this argument 'time'.")

        return positional_args, keyword_args

    def _compute_value_args(self):
        """Compute the lists of arguments."""
        self.positional_args, self.keyword_args = self.compute_args(self.value)

    @property
    def value(self):
        """The value."""
        return self._value

    @value.setter
    def value(self, value):
        self._value = value
        self._compute_value_args()

    @property
    def inverse_value(self):
        """The inverse value"""
        if callable(self.value):
            parameters = signature(self.value).parameters
            defaults = [parameters.get(kwarg) for kwarg in self.keyword_args]
            if len(self.positional_args) == 0:
                variables_head = ', '.join(f'{v}={d.default}' for v, d in zip(self.keyword_args, defaults))
                variables_body = ', '.join(f'{v}={v}' for v in self.keyword_args)
            else:
                variables_head = ', '.join([', '.join(self.positional_args),
                                            ', '.join(f'{v}={d.default}' for v, d in
                                                      zip(self.keyword_args, defaults))])
                variables_body = ', '.join([', '.join(self.positional_args),
                                            ', '.join(f'{v}={v}' for v in self.keyword_args)])
            head = f'def inverse_fun({variables_head}):'
            body = f'return 1/value({variables_body})'
            local_dict = {}
            exec(f'{head}{body}', {'value': self._value}, local_dict)  # pylint: disable=exec-used

            return local_dict['inverse_fun']
        return 1 / self.value

    @property
    def is_homogeneous(self) -> bool:
        """
        Returns True if the value is homogeneous and False if not

        Returns
        -------
        bool
            True: value is homogeneous. False: value is not homogeneous.

        """
        if not callable(self.value):
            return True
        return len(self.positional_args) == 0

    @property
    def is_linear(self) -> bool:
        """
        Returns True if the value is linear and False if not

        Returns
        -------
        bool
            True: value is linear. False: value is not linear.

        """
        if not callable(self.value):
            return True
        return len(set(self.keyword_args) - {'time'}) == 0

    @property
    def is_time_dependent(self) -> bool:
        """
        Returns True if the value is time dependent and False if not

        Returns
        -------
        bool
            True: value is time dependent. False: value is not time dependent.

        """
        return 'time' in self.keyword_args

    @property
    def is_constant(self) -> bool:
        """Returns True if the value is constant, i.e. linear, homogeneous and not time dependent

        Returns
        -------
        constant : bool
        """
        return all((self.is_homogeneous, self.is_linear, not self.is_time_dependent))

    def update(self, name: str, val: Any):
        """Updates the arguments of the value, if possible.

        Parameters
        ----------
        name : str
            The name of the argument of value to update
        val : Any
            The new value
        """
        if name not in self.keyword_args:
            return
        self._value.__setattr__(name, val)
