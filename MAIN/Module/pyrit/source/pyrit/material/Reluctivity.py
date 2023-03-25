# coding=utf-8
"""Implementation of Reluctivity.

.. sectionauthor:: Bundschuh
"""

from . import MatProperty


class Reluctivity(MatProperty):
    r"""Class that represents the Reluctivity :math:`\nu=\mu^{-1}`"""

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
