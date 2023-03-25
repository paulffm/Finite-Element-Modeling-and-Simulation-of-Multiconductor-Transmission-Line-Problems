# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 13:01:27 2021

.. sectionauthor:: leissner
"""

import numpy as np
import gmsh

from . import PrimitiveGeo, Point


# from pyrit.geometry.Extrude import Extrude
# from pyrit.geometry.Pipes import Pipes
# from pyrit.geometry.Boolean import Boolean


class Curve(PrimitiveGeo):
    """A class representing a curve."""

    __slots__ = ('_points', '_num_points', '_reverse', '_tag')

    def __init__(self, *points: Point, physical_tag: int = None):
        """
        Constructor of the class Curve

        Parameters
        ----------
        *points : Point
            A set of points that are part of the curve
        physical_tag : int, optional
            The physical tag of the curve. The default is None.

        Returns
        -------
        None.

        """
        super().__init__(1, physical_tag)
        self._points = points
        self._num_points = len(points)
        self._reverse = False
        self._tag = None

    @property
    def dim(self) -> int:
        """
        Get the dimension of the Curve.

        Returns
        -------
        int
            The dimension

        """
        return self._dim

    @property
    def points(self):
        """
        Get the points of a curve.

        Returns
        -------
        Tuple
            Tuple in which each element represents a point.

        """
        return self._points

    @property
    def num_points(self) -> int:
        """
        Get the number of points that are part of the Curve.

        Returns
        -------
        int
            Number of points part of the curve.

        """
        return self._num_points

    @property
    def added(self) -> bool:
        """
        Give True if the entity has beed added to gmsh and False if not.

        Returns
        -------
        bool
            True if added in gmsh False if not

        """
        return self.tag is not None

    @property
    def reversed(self):
        """After using this property, next time getting the (existing) tag, the tag is negative."""
        self._reverse = True
        return self

    @property
    def tag(self):
        """Get the tag."""
        if self._reverse:
            if self._tag is not None:
                self._reverse = False
                return -1 * self._tag
        return self._tag

    @tag.setter
    def tag(self, tag):
        """Set the tag."""
        self._tag = tag

    def reset_tag(self) -> None:
        self.tag = None
        for point in self._points:
            point.tag = None

    def remove(self, recursive=True):
        """
        Remove the Curve from GMSH.

        Parameters
        ----------
        recursive : bool, optional
            If True, all entities that define this entity are removed as well. Default is True

        Returns
        -------
        None.

        """
        if self.added:
            gmsh.model.occ.remove([(self._dim, self.tag)], recursive=recursive)
            if recursive:
                self.reset_tag()
            else:
                self.tag = None
        else:
            raise Warning("This point has not been added and can thus not be removed")

    def translate(self, dx: float, dy: float, dz: float):
        """
        Translate the Curve.

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
        # translating the points of the curve
        for i in self._points:
            i.translate(dx, dy, dz)

        # translating the curve in gmsh
        if self.added:
            gmsh.model.occ.translate([(self._dim, self.tag)],
                                     dx, dy, dz)

    def rotate(self, angle_x: float, angle_y: float, angle_z: float,
               rotation_origin: np.ndarray = None):
        """
        Rotate the Curve.

        Note the order of the rotation: first rotate around x-axis,
        then rotate around y-axis, the rotate around z-axis.

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
        # rotating the points of the curve
        for i in self._points:
            i.rotate(angle_x, angle_y, angle_z, rotation_origin)

        # rotating the curve in gmsh
        if self.added:
            if rotation_origin[0] is not None:
                x = rotation_origin[0]
            else:
                x = 0
            if rotation_origin[1] is not None:
                y = rotation_origin[1]
            else:
                y = 0
            if rotation_origin[2] is not None:
                z = rotation_origin[2]
            else:
                z = 0
            gmsh.model.occ.rotate([(self._dim, self.tag)],
                                  x, y, z, 1, 0, 0, angle_x)
            gmsh.model.occ.rotate([(self._dim, self.tag)],
                                  x, y, z, 0, 1, 0, angle_y)
            gmsh.model.occ.rotate([(self._dim, self.tag)],
                                  x, y, z, 0, 0, 1, angle_z)

    def add_to_gmsh(self):
        raise Warning("Curve by itself cannot be added in GMSH")
