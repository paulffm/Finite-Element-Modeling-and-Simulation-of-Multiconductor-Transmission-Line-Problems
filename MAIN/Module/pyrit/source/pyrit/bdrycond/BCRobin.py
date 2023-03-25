# -*- coding: utf-8 -*-
"""
Created on Thu Mai 27 2021

Contains the class BCRobin
.. sectionauthor:: menzenbach
"""

from typing import Callable, Any, Union
from numpy import ndarray

from . import BCUnary


class BCRobin(BCUnary):
    r"""
    Class of BCRobin with super class BCUnary.

    Notes
    -----
    For a scalar potential, let :math:`\alpha` be the coefficient of the Dirichlet part, :math:`\beta` the
    coefficient of the Neumann part and :math:`g` the value. Then it is :math:`\alpha u + \beta \epsilon\frac{\partial
    u}{\partial n} = g` the Robin boundary condition.
    """

    def __init__(self, coefficient_dirichlet: Union[Callable[[Any], float], ndarray, float],
                 coefficient_neumann: Union[Callable[[Any], float], ndarray, float],
                 value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = '') -> None:
        """
        Constructor of class BCRobin.

        Parameters
        ----------
        coefficient_dirichlet: Union[Callable[[Any], float], ndarray, float]
            Coefficient of Dirichlet part
        coefficient_neumann: Union[Callable[[Any], float], ndarray, float]
            Coefficient of Neumann part
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents right hand side of the boundary condition
        name : str, optional
            Name of the BC, by default ''
        """
        self.positional_args_dirichlet = []
        self.keyword_args_dirichlet = []
        self.coefficient_dirichlet = coefficient_dirichlet
        self.coefficient_neumann = coefficient_neumann
        super().__init__(value, name=name)

    @property
    def coefficient_dirichlet(self):
        """The value."""
        return self._coefficient_dirichlet

    @coefficient_dirichlet.setter
    def coefficient_dirichlet(self, coefficient_dirichlet):
        self._coefficient_dirichlet = coefficient_dirichlet
        self.positional_args_dirichlet, self.keyword_args_dirichlet = self.compute_args(coefficient_dirichlet)

    @property
    def coefficient_neumann(self):
        """The value."""
        return self._coefficient_neumann

    @coefficient_neumann.setter
    def coefficient_neumann(self, coefficient_neumann):
        self._coefficient_neumann = coefficient_neumann
        self.positional_args_neumann, self.keyword_args_neumann = self.compute_args(coefficient_neumann)

    def _is_homogeneous(self, name_call: str, name_args: str) -> bool:
        """General function for checking for homogeneity.

        Parameters
        ----------
        name_call : str
            Name of the argument to check.
        name_args : str
            Name of thr args-list.

        Returns
        -------
        is_homogeneous : bool
            Weather the argument of homogeneous or not
        """
        if not callable(self.__getattribute__(name_call)):
            return True
        return len(self.__getattribute__(name_args)) == 0

    def _is_linear(self, name_call: str, name_args: str) -> bool:
        """General function for checking for linearity.

        Parameters
        ----------
        name_call : str
            Name of the argument to check.
        name_args : str
            Name of thr args-list.

        Returns
        -------
        is_linear : bool
            Weather the argument is linear or not
        """
        if not callable(self.__getattribute__(name_call)):
            return True
        return len(set(self.__getattribute__(name_args)) - {'time'}) == 0

    def _is_time_dependent(self, name: str) -> bool:
        """General function for checking for time dependency.

        Parameters
        ----------
        name : str
            Name of the keyword args argument.

        Returns
        -------
        is_time_dependent : bool
            Weather the argument is time dependent or not.
        """
        return 'time' in self.__getattribute__('keyword_args' + name)

    @property
    def is_homogenous_dirichlet(self) -> bool:
        """
        Returns True if the dirichlet value is homogeneous and False if not

        Returns
        -------
        bool
            True: dirichlet value is homogeneous. False: dirichlet value is not homogeneous.

        """
        return self._is_homogeneous('_coefficient_dirichlet', 'positional_args_dirichlet')

    @property
    def is_homogeneous_neumann(self) -> bool:
        """
        Returns True if the neumann value is homogeneous and False if not

        Returns
        -------
        bool
            True: neumann value is homogeneous. False: neumann value is not homogeneous.

        """
        return self._is_homogeneous('_coefficient_neumann', 'positional_args_neumann')

    @property
    def is_homogeneous_value(self) -> bool:
        """
        Returns True if the value is homogeneous and False if not

        Returns
        -------
        bool
            True: value is homogeneous. False: value is not homogeneous.

        """
        return self._is_homogeneous('_value', 'positional_args')

    @property
    def is_homogeneous(self) -> bool:
        """
        Returns True if the value, dirichlet and neumann value are homogeneous and False if not

        Returns
        -------
        bool
            True: all are homogeneous. False: not all are homogeneous.

        """
        return all((self.is_homogenous_dirichlet, self.is_homogeneous_neumann, self.is_homogeneous_value))

    @property
    def is_linear_dirichlet(self) -> bool:
        """
        Returns True if the dirichelt value is linear and False if not

        Returns
        -------
        bool
            True: dirichlet value is linear. False: dirichlet value is not linear.

        """
        return self._is_linear('_coefficient_dirichlet', 'keyword_args_dirichlet')

    @property
    def is_linear_neumann(self) -> bool:
        """
        Returns True if the neumann value is linear and False if not

        Returns
        -------
        bool
            True: neumann value is linear. False: neumann value is not linear.

        """
        return self._is_linear('_coefficient_neumann', 'keyword_args_neumann')

    @property
    def is_linear_value(self) -> bool:
        """
        Returns True if the value is linear and False if not

        Returns
        -------
        bool
            True: value is linear. False: value is not linear.

        """
        return self._is_linear('value', 'keyword_args')

    @property
    def is_linear(self) -> bool:
        """
        Returns True if the value, dirichlet and neumann value is linear and False if not

        Returns
        -------
        bool
            True: all are linear. False: not all are not linear.

        """
        return all((self.is_linear_dirichlet, self.is_linear_neumann, self.is_linear_value))

    @property
    def is_time_dependent_dirichlet(self) -> bool:
        """
        Returns True if the dirichlet value is time dependent and False if not

        Returns
        -------
        bool
            True: dirichlet value is time dependent. False: dirichlet value is not time dependent.
        """
        return self._is_time_dependent('_dirichlet')

    @property
    def is_time_dependent_neumann(self) -> bool:
        """
        Returns True if the neumann value is time dependent and False if not

        Returns
        -------
        bool
            True: neumann value is time dependent. False: neumann value is not time dependent.
        """
        return self._is_time_dependent('_neumann')

    @property
    def is_time_dependent_value(self) -> bool:
        """
        Returns True if the value is time dependent and False if not

        Returns
        -------
        bool
            True: value is time dependent. False: value is not time dependent.
        """
        return self._is_time_dependent('')

    @property
    def is_time_dependent(self) -> bool:
        """
        Returns True if the value, dirichlet and neumann value is time dependent and False if not

        Returns
        -------
        bool
            True: all are time dependent. False: not all are time dependent.
        """
        return all((self.is_time_dependent_dirichlet, self.is_time_dependent_neumann, self.is_time_dependent_value))

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
        for kw, arg in zip((self.keyword_args, self.keyword_args_dirichlet, self.keyword_args_neumann),
                           ('_value', '_coefficient_dirichlet', '_coefficient_neumann')):
            if name not in kw:
                continue
            self.__getattribute__(arg).__setattr__(name, val)
