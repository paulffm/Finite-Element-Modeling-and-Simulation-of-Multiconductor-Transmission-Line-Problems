# -*- coding: utf-8 -*-
"""
Created on Thu Mai 27 2021

Contains the class BCDirichlet

.. sectionauthor:: menzenbach
"""

from typing import Callable, Any, Union
from numpy import ndarray

from . import BCUnary


class BCDirichlet(BCUnary):
    r"""Class of BCDirichlet with super class BCUnary.

    Notes
    -----
    A Dirichlet boundary condition defines a fixed value :math:`\gamma` on the boundary :math:`\Gamma_{\mathrm{bnd}}`,
    e.g. a fixed potential in V or a temperature in K, yielding: :math:`\phi|_{\Gamma_{\mathrm{bnd}}}=\gamma`.
    """

    def __init__(self, value: Union[Callable[[Any], Union[ndarray, float]], ndarray, float], name: str = '') -> None:
        """
        Constructor of class BCDirichlet.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the BC
        name : str, optional
            Name of the BC, by default ''
        """
        super().__init__(value, name=name)
