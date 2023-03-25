# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:18:22 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty


class ThermalConductivity(MatProperty):
    r"""Class that represents the thermal conductivity :math:`\lambda`"""

    def __init__(self, value, name: str = ''):
        """
        Class that represents the thermal conductivity

        Parameters
        ----------
        value :
            Value of the ThermalConductivity (see source.material.MatProperty).
        name : str, optional
            Name of the ThermalConductivity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
