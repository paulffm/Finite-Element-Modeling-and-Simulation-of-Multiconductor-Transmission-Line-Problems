# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:17:45 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty
from . import Conductivity


class DifferentialConductivity(MatProperty):
    r"""Class that represents the differential conductivity :math:`\boldsymbol{\frac{\mathrm{d} J}{\mathrm{d} E}}`.

    In case of a field-dependent conductivity :math:`\sigma(\boldsymbol{E})`, it is defined as

    .. math::
        \boldsymbol{\frac{\mathrm{d} J}{\mathrm{d} E}} = \sigma(\boldsymbol{E}) \cdot \mathrm{unit\_tensor} + 2 \cdot
        \frac{\mathrm{d} \sigma}{\mathrm{d} E^2}(\boldsymbol{E}) \cdot \boldsymbol{E}\boldsymbol{E}^{\mathrm{T}}.

    The definition is supported by :py:func:`pyrit.toolbox.MaterialLibrary.get_differential_tensor`, an example
    is given in tutorials/nonlinear_resistor.py.
    """

    basis_class = Conductivity

    def __init__(self, value, name: str = ''):
        """
        Constructor of DifferentialConductivity

        Parameters
        ----------
        value :
            Value of the differential conductivity (see source.material.MatProperty).
        name : str, optional
            Name of the differential conductivity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
