# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 12:17:45 2021

.. sectionauthor:: bundschuh
"""

from . import MatProperty, Permittivity


class DifferentialPermittivity(MatProperty):
    r"""Class that represents the differential permittivity :math:`\boldsymbol{\frac{\mathrm{d} D}{\mathrm{d} E}}`.

    In case of a field-dependent permittivity :math:`\varepsilon(\boldsymbol{E})`, it is defined as

    .. math::
        \boldsymbol{\frac{\mathrm{d} D}{\mathrm{d} E}} = \varepsilon(\boldsymbol{E}) \cdot \mathrm{unit\_tensor}
        + 2 \cdot \frac{\mathrm{d} \varepsilon}{\mathrm{d} E^2}(\boldsymbol{E}) \cdot
        \boldsymbol{E}\boldsymbol{E}^{\mathrm{T}}.

    The definition is supported by :py:func:`pyrit.toolbox.MaterialLibrary.get_differential_tensor`, an example
    is given in tutorials/nonlinear_resistor.py.
    """

    basis_class = Permittivity

    def __init__(self, value, name: str = ''):
        """
        Constructor of differential permittivity

        Parameters
        ----------
        value :
            Value of the differential permittivity (see source.material.MatProperty).
        name : str, optional
            Name of the differential permittivity. The default is ''.

        Returns
        -------
        None.

        """
        super().__init__(value, name)
