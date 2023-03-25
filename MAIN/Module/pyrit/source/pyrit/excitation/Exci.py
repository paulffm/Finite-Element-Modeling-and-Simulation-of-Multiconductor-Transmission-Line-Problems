# coding=utf-8
"""
Implementation of Exci.

.. sectionauthor:: Bundschuh
"""
# Linting is here disables because of an error in pylint https://github.com/PyCQA/pylint/issues/214 with duplicate code.
# The file was linted.
# pylint: disable-all

from abc import ABC, abstractmethod
from inspect import signature
from typing import Any, Callable, Union
from numpy import ndarray

from pyrit import ValueHandler


class Exci(ValueHandler, ABC):
    """Abstract class for excitations."""

    __ID_count = 0  # Provides unique ID for the Mat

    @abstractmethod
    def __init__(self, value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = ''):
        """
        Constructor of Exci

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the Exci
        name : str, optional
            A name of the Exci
        """
        super().__init__(value)
        self._ID = Exci.__ID_count
        self.name = name

        # Increasing the id count
        Exci.__ID_count += 1

    @property
    def ID(self) -> int:
        """
        The unique ID of the Exci.

        Returns
        -------
        ID: int
        """
        return self._ID
