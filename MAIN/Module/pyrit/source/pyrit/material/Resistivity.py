# coding=utf-8
"""Implementation of Resistivity

.. sectionauthor:: Bundschuh
"""
from pyrit.material.MatProperty import MatProperty


class Resistivity(MatProperty):
    r"""Class representing the resistivity :math:`\rho=\sigma^{-1}`"""

    def __init__(self, value, name: str = ''):
        """
        Constructor of Resistivity

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
