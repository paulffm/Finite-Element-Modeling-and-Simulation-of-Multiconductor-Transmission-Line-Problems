# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 2021

Contains class Circle

.. sectionauthor:: leissner
"""

import gmsh
import numpy as np


from . import Point, Curve


class Circle(Curve):
    """Class representing a circle."""

    __slots__ = ("center_x", "center_y", "center_z", "__radius",
                 "angle1", "angle2", "angle_x", "angle_y", "angle_z")

    def __init__(self, center_x: float, center_y: float, center_z: float,
                 radius: float, angle1: float = 0, angle2: float = 2 * np.pi,
                 angle_x: float = 0, angle_y: float = 0,
                 angle_z: float = 0, physical_tag: int = None):
        """
        Create an object of the class Circle.

        Parameters
        ----------
        center_x : float
            x-coordinate of the center point.
        center_y : float
            y-coordinate of the center point.
        center_z : float
            z-coordinate of the center point.
        radius : float
            Radius of the circle.
        angle1 : float, optional
            Angle at which the circle arc begins. The default is None.
        angle2 : float, optional
            Angle at which the circle arc ends. The default is None.
        angle_x : float, optional
            Angle to rotate the circle around it's x-axis. The default is None.
        angle_y : float, optional
            Angle to rotate the circle around it's y-axis. The default is None.
        angle_z : float, optional
            Angle to rotate the circle around it's z-axis. The default is None.
        physical_tag : int, optional
            Physical_tag of the circle. The default is None.

        Returns
        -------
        None.

        """
        self.center_x = center_x
        self.center_y = center_y
        self.center_z = center_z
        self.radius = radius
        self.angle1 = angle1
        self.angle2 = angle2
        self.angle_x = angle_x
        self.angle_y = angle_y
        self.angle_z = angle_z

        # start_point and end_point are only created in the xy-plane by gmsh
        # values of these points are later assigned when added in gmsh
        # end_point only exists if it is not a closed circle
        if np.isclose(angle1 + 2 * np.pi, angle2):
            start_point = Point(None, None, None)
            super().__init__(start_point, physical_tag=physical_tag)
        else:
            start_point = Point(None, None, None)
            end_point = Point(None, None, None)
            super().__init__(start_point, end_point, physical_tag=physical_tag)

    @property
    def start_point(self) -> Point:
        """
        Get the start point of the circle. At this point the circle arc begins.

        Returns
        -------
        Point
            Start point of the circle.

        """
        return self.points[-1]

    @property
    def end_point(self) -> Point:
        """
        Get the end_point of the circle. At this point the circle arc ends.

        Returns
        -------
        Point
            End point of the circle.

        """
        return self.points[0]

    @property
    def radius(self) -> float:
        """
        Get the radius of the circle.

        Returns
        -------
        float
            Radius of the circle.

        """
        return self.__radius

    @radius.setter
    def radius(self, radius: float):
        """
        Set the radius of the circle.

        Parameters
        ----------
        radius : float
            Value of the new radius.

        Raises
        ------
        TypeError
            If radius is not int or float and if radius is less or equal to 0.
        """
        if isinstance(radius, float) or isinstance(radius, int) and radius > 0:
            self.__radius = radius
        else:
            raise TypeError(
                "Unable to set radius. Radius must be a float greater than 0.")

    def add_to_gmsh(self):
        """
        Add the circle in gmsh.

        Raises
        ------
        Warning
            If the circle has been already added in GMSH.

        Returns
        -------
        None.

        """
        if self.added:  # pylint: disable=no-else-raise
            raise Warning("The circle has already been added in GMSH!")
        else:
            # Add the circle on the xy-plane with the center point
            # (self.center_x, self.center_y, self.center_z):
            self.tag = gmsh.model.occ.addCircle(self.center_x, self.center_y,
                                                self.center_z, self.radius,
                                                angle1=self.angle1,
                                                angle2=self.angle2)

            gmsh.model.occ.synchronize()  # run CAD-kernel

            # Add tag of a point from gmsh to each point:
            # gmsh returns the endpoint first and then the startpoint
            entities_0D = gmsh.model.occ.getEntities(0)

            start_point_coord = gmsh.model.getValue(0, entities_0D[len(entities_0D) - 1][1], [])
            self.points[0].x = start_point_coord[0]
            self.points[0].y = start_point_coord[1]
            self.points[0].z = start_point_coord[2]
            self.points[0].tag = entities_0D[len(entities_0D) - 1][1]

            if self.num_points == 2:
                self.points[1].tag = entities_0D[len(entities_0D) - 2][1]

                end_point_coord = gmsh.model.getValue(0, entities_0D[len(entities_0D) - 2][1], [])
                self.points[1].x = end_point_coord[0]
                self.points[1].y = end_point_coord[1]
                self.points[1].z = end_point_coord[2]

            # 3D placement of the circle:
            super().rotate(self.angle_x, self.angle_y, self.angle_z,
                           np.array([self.center_x, self.center_y,
                                     self.center_z]))
