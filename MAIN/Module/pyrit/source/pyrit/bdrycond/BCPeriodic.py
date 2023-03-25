# -*- coding: utf-8 -*-
"""
Created on Tue June 15 2021

Contains the class BCPeriodic
.. sectionauthor:: menzenbach
"""


from . import BCBinary


class BCPeriodic(BCBinary):  # todo: Add physical meaning and possible usages.
    """Class of BCPeriodic with super class BCBinary.

    The value is set to 1.
    """

    def __init__(self, name: str = ''):
        """
        Constructor of class BCPeriodic.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the BC
        name : str, optional
            Name of the BC, by default ''
        """
        super().__init__(1, name=name)
