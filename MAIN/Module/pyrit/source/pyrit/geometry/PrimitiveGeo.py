# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 10:16:23 2021

Contains the abstract class PrimitiveGeo

.. sectionauthor:: bundschuh
"""

from abc import ABC, abstractmethod
import numpy as np


class PrimitiveGeo(ABC):
    """Abstract class for primitive geometry object that can be added directly to GMSH."""

    __slots__ = ('tag', '_dim', 'physical_tag')

    @abstractmethod
    def __init__(self, dim: int, physical_tag: int, tag=None):
        self._dim = dim
        if not isinstance(physical_tag, int) and physical_tag is not None:
            raise TypeError("physical_tag needs to be an integer or None!")
        self.physical_tag = physical_tag
        self.tag = None

    @property
    @abstractmethod
    def dim(self) -> int:
        """
        Get the dimension of the PrimitiveGeo object

        The following dimensions exist:
            0: Point
            1: Line
            2: Surface
            3: Volume

        Returns
        -------
        int
            The dimension

        """

    @property
    @abstractmethod
    def added(self) -> bool:
        """
        Gives True if the entity has beed added to gmsh and False if not

        Returns
        -------
        bool
            True if added in gmsh False if not

        """

    @abstractmethod
    def reset_tag(self) -> None:
        """
        Sets its own tag and the tags of all contained entities to None.

        Returns
        -------
        None
        """

    @abstractmethod
    def add_to_gmsh(self):
        """
        Add the PrimitiveGeo object in GMSH

        Returns
        -------
        None.

        """

    @abstractmethod
    def remove(self, recursive=True):
        """
        Remove the PrimitiveGeo object from GMSH

        Parameters
        ----------
        recursive : bool, optional
            If True, all entities that define this entity are removed as well. Default is True

        Returns
        -------
        None.

        """

    @abstractmethod
    def translate(self, dx: float, dy: float, dz: float):
        """
        Translate the PrimitiveGeo object.

        Parameters
        ----------
        dx : float
            Translation in x-direction
        dy : float
            Translation in y-direction
        dz : float
            Translation in z-direction

        Returns
        -------
        None.

        """

    @abstractmethod
    def rotate(self, angle_x: float, angle_y: float, angle_z: float,
               rotation_origin: np.ndarray = None):
        """
        Rotate the PrimitiveGeo object.

        Parameters
        ----------
        angle_x : float
            Rotation angle around x-axis in radian
        angle_y : float
            Rotation angle around y-axis in radian
        angle_z : float
            Rotation angle around z-axis in radian
        rotation_origin : np.ndarray, optional
            The origin around which is rotated. By default (0,0,0) is taken

        Returns
        -------
        None.

        """
