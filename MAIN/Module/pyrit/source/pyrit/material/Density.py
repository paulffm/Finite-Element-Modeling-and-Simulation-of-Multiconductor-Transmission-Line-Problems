# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:19:25 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty


class Density(MatProperty):
    """Class that represents the density"""

    def __init__(self, value, name: str = ''):
        """
        Constructor of Density

        Parameters
        ----------
        value :
            Value of the Density (see source.material.MatProperty).
        name : str, optional
            Name of the Density. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
