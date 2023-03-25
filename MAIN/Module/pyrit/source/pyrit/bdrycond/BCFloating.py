# -*- coding: utf-8 -*-
"""
Created on Thu Mai 27 2021

Contains the class BCFloating
.. sectionauthor:: menzenbach
"""


from . import BCUnary


class BCFloating(BCUnary):
    r"""Class of BCFloating with super class BCUnary.

    Notes
    -----
    The value of the potential on this boundary is constant, but unknown beforehand.
    Thus, a fixed *unknown* value :math:`\gamma` is defined on the boundary :math:`\Gamma_{\mathrm{bnd}}`,
    e.g. a fixed potential in V or a temperature in K, yielding :math:`\phi|_{\Gamma_{\mathrm{bnd}}}=\gamma`.
    This expression is inserted in the system of equations and solved accordingly.
    """

    def __init__(self, name: str = '') -> None:
        """
        Constructor of class BCFloating.

        Parameters
        ----------
        value : Union[Callable[[Any], Union[ndarray, float]], ndarray, float]
            A function or value that represents the BC
        name : str, optional
            Name of the BC, by default ''
        """
        super().__init__(None, name=name)
