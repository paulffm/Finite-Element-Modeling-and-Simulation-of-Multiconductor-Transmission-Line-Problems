# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:17:45 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty


class Conductivity(MatProperty):
    r"""Class that represents the conductivity :math:`\sigma` in S/m."""

    def __init__(self, value, name: str = ''):
        """
        Constructor of Conductivity

        Parameters
        ----------
        value :
            Value of the Conductivity (see source.material.MatProperty).
        name : str, optional
            Name of the Conductivity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
