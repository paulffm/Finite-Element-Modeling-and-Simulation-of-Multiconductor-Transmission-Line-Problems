# -*- coding: utf-8 -*-
"""
Created on Tue June 15 2021

Contains the class BCAntiPeriodic
.. sectionauthor:: menzenbach
"""


from . import BCBinary


class BCAntiPeriodic(BCBinary):
    """Class of BCAntiPeriodic with super class BCBinary.

    The value is set to -1. #todo: Add physical meaning and
    possible usages.
    """

    def __init__(self, name: str = ''):
        """
        Constructor of class BCAntiPeriodic.

        Parameters
        ----------
        name : str, optional
            Name of the BC, by default ''
        """
        super().__init__(-1, name=name)
