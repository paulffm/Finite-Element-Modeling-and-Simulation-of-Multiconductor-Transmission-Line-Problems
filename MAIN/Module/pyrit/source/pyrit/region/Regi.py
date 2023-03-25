# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 2021

Contains the class Regi

.. sectionauthor:: menzenbach
"""

from typing import Optional


class Regi:
    """
    Class for Regi.

    Attributes
    ----------
    id : int
        Unique ID for Regi
    mat : int
        Material ID
    bc: int
        Bounding Condition ID
    exci: int
        Excitation ID
    """

    def __init__(self, ID: int, dim: int, mat: Optional[int] = None, bc: Optional[int] = None,
                 exci: Optional[int] = None):
        """
        Constructor of class Regi

        Parameters
        ----------
        ID : int
            Unique ID for Regi
        dim: int
            Dimension of the Regi
        mat : int, optional
            Material ID, by default None
        bc : int, optional
            Bounding Condition ID, by default None
        exci : int, optional
            Excitation ID, by default None
        """
        self.ID = ID
        self.dim = dim
        self.mat = mat
        self.bc = bc
        self.exci = exci

    def __eq__(self, other):
        # pylint: disable=line-too-long
        return self.ID == other.ID and self.dim == other.dim and self.mat == other.mat and self.bc == other.bc and self.exci == other.exci
