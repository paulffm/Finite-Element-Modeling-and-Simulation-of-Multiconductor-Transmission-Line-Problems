# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 2021

Contains the abstract class BCUnary

.. sectionauthor:: menzenbach
"""

from abc import ABC, abstractmethod
from typing import Callable, Any, Union
from numpy import ndarray

from . import BC


class BCUnary(BC, ABC):  # todo: Add physical meaning and
    """
    Abstract class of a unary boundary condition object.

    possible usages.

    Attributes
    ----------
    id : int
        Unique identifier for each object
    name : str
        A name for the boundary condition
    regi : int
        ID of the region in which the boundary is created
    idxoffset : double
        Offset regarding the degrees of freedom in the solution vector
    data : complex double vector
        Boundary data (constant part)
    timefun : function handle
        Boundary data (part depending on time)
    currentval: double vector
        Current value of the boundary data
    """

    @abstractmethod
    def __init__(self, value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = '') -> None:
        """
        Constructor of abstract class BCUnary.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            The value of the BC. A function or a constant
        name : str
            A name for the boundary condition
        """
        super().__init__(value, name)
