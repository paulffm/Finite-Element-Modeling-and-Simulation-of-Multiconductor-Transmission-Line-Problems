# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 11:03:01 2021

Contains the class Point

.. sectionauthor:: bundschuh
"""

import numpy as np
import gmsh

from . import PrimitiveGeo


class Point(PrimitiveGeo):
    """A class representing a point."""

    __slots__ = ('x', 'y', 'z', '__mesh_size')

    def __init__(self, x: float, y: float, z: float, mesh_size: float = None, physical_tag: int = None):
        """
        Constructor of the class Point

        Parameters
        ----------
        x : float
            x-coordinate of the point
        y : float
            y-coordinate of the point
        z : float
            z-coordinate of the point
        mesh_size : float, optional
            A point specific mesh size factor. The default is None.
        physical_tag: int, optional
            Saves the tag of the physical group

        Returns
        -------
        None.

        """
        super().__init__(0, physical_tag)
        self.x = x
        self.y = y
        self.z = z
        self.mesh_size = mesh_size

    def __repr__(self):
        return f"Point with tag '{self.tag}' at coordinates ({self.x},{self.y},{self.z})"

    @property
    def dim(self) -> int:
        return self._dim

    @property
    def mesh_size(self) -> float:
        """Returns the mesh size at this point."""
        return self.__mesh_size

    @property
    def added(self) -> bool:
        return self.tag is not None

    @mesh_size.setter
    def mesh_size(self, mesh_size: float):
        """
        Sets the mesh size

        Parameters
        ----------
        mesh_size : float
            The new mesh size

        Raises
        ------
        ValueError
            When the new mesh size is not a positive number

        Returns
        -------
        None.

        """
        if mesh_size is None:
            self.__mesh_size = mesh_size
        elif mesh_size > 0:
            self.__mesh_size = mesh_size
        else:
            raise ValueError(f"The mesh_size needs to be a positive number. It was {mesh_size}")

    def reset_tag(self) -> None:
        self.tag = None

    def add_to_gmsh(self):
        if not self.added:
            if self.__mesh_size is not None:
                tag = gmsh.model.occ.addPoint(self.x, self.y, self.z, self.__mesh_size)
            else:
                tag = gmsh.model.occ.addPoint(self.x, self.y, self.z)
            self.tag = tag
        else:
            raise Warning("This Point has already been adden in gmsh. It will not be added again.")

    def remove(self, recursive=True):
        if self.added:
            gmsh.model.occ.remove([(self._dim, self.tag)])
            gmsh.model.removeEntities([(self._dim, self.tag)])
            self.tag = None
        else:
            raise Warning("This point has not been added and can thus not be removed")

    def translate(self, dx: float, dy: float, dz: float):
        self.x += float(dx)
        self.y += float(dy)
        self.z += float(dz)

    def rotate(self, angle_x: float, angle_y: float, angle_z: float,
               rotation_origin: np.ndarray = None):
        if rotation_origin is None:
            rotation_origin = np.array([0, 0, 0])

        p = np.array([self.x, self.y, self.z])
        p_prime = p - rotation_origin  # Transform

        # Create rotation matrices for each axis
        rotation_x = np.array(
            [[1, 0, 0], [0, np.cos(angle_x), -1 * np.sin(angle_x)], [0, np.sin(angle_x), np.cos(angle_x)]])
        rotation_y = np.array(
            [[np.cos(angle_y), 0, np.sin(angle_y)], [0, 1, 0], [-1 * np.sin(angle_y), 0, np.cos(angle_y)]])
        rotation_z = np.array(
            [[np.cos(angle_z), -1 * np.sin(angle_z), 0], [np.sin(angle_z), np.cos(angle_z), 0], [0, 0, 1]])

        rotation = rotation_z @ rotation_y @ rotation_x  # Compute the overall rotation matrix

        p_rotated_prime = rotation @ p_prime  # Rotate

        p_rotated = p_rotated_prime + rotation_origin  # Transform back

        (self.x, self.y, self.z) = p_rotated.tolist()  # Update coordinates
