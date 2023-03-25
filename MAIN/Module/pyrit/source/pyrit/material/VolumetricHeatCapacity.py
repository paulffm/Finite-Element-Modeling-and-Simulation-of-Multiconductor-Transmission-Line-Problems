# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:19:02 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty


class VolumetricHeatCapacity(MatProperty):
    r"""Class that represents the volumetric heat capacity :math:`c_v = \rho c_p` (German: Wärmespeicherzahl)"""

    def __init__(self, value, name: str = ''):
        """
        Constructor of VolumetricHeatCapacity

        Parameters
        ----------
        value :
            Value of the VolumetricHeatCapacity (see source.material.MatProperty).
        name : str, optional
            Name of the VolumetricHeatCapacity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
