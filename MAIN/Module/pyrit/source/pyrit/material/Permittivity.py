# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:02:36 2021

.. sectionauthor:: bundschuh
"""


from . import MatProperty


class Permittivity(MatProperty):
    r"""Class that represents the Permittivity :math:`\varepsilon=\varepsilon_{\mathrm{r}} \varepsilon_0`"""

    def __init__(self, value, name: str = ''):
        """
        Constructor of Permittivity

        Parameters
        ----------
        value :
            Value of the Permittivity (see source.material.MatProperty).
        name : str, optional
            Name of the Permittivity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
