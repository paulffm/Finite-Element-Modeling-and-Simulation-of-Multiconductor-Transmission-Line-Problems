# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 2021

Contains the abstract class BC

.. sectionauthor:: menzenbach
"""
# Linting is here disables because of an error in pylint https://github.com/PyCQA/pylint/issues/214 with duplicate code.
# The file was linted.
# pylint: disable-all
from abc import ABC, abstractmethod
from inspect import signature
from typing import Callable, Any, Union
from numpy import ndarray

from pyrit import ValueHandler


class BC(ValueHandler, ABC):
    """Abstract class of a boundary condition object."""

    __ID_count = 0      # Provides an unique ID for the BC

    @abstractmethod
    def __init__(self, value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = ''):
        """
        Constructor of the abstract class BC.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the BC
        name : str, optional
            A name of the BC
        """
        super().__init__(value)
        self._ID = BC.__ID_count
        self.name = name

        # Increasing the id count
        BC.__ID_count += 1

    @property
    def ID(self) -> int:
        """
        The unique ID of the BC

        Returns
        -------
        ID: int
        """
        return self._ID
