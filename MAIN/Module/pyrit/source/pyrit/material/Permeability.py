# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:17:11 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty


class Permeability(MatProperty):
    r"""Class that represents the permeability :math:`\mu=\mu_{\mathrm{r}} \mu_0`"""

    def __init__(self, value, name: str = ''):
        """
        Constructor of Permeability

        Parameters
        ----------
        value :
            Value of the Permeability (see source.material.MatProperty).
        name : str, optional
            Name of the Permeability. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
