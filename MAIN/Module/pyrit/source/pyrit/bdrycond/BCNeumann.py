# -*- coding: utf-8 -*-
"""
Created on Thu Mai 27 2021

Contains the class BCNeumann
.. sectionauthor:: menzenbach
"""

from typing import Callable, Any, Union
from numpy import ndarray

from . import BCUnary


class BCNeumann(BCUnary):
    r"""Class of BCNeumann with super class BCUnary.

    Notes
    -----
    The Neumann boundary condition specifies the value :math:`g` of the derivative of the unknown :math:`u` with
    respect to the normal direction, applied at the boundary of the domain, :math:`\frac{\partial u}{\partial \vec{
    n}} = g`. It is the natural (or implicitly defined) boundary condition in a FE setting.
    """

    def __init__(self, value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = '') -> None:
        """
        Constructor of class BCNeumann.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the BC
        name : str, optional
            Name of the BC, by default ''
        """
        super().__init__(value, name=name)
