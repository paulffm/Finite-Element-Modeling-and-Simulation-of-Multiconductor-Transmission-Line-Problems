# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 2021

Contains the class Ellipse

.. sectionauthor:: menzenbach
"""

import warnings
import math
import gmsh
import numpy as np


from . import Curve, Point


class Ellipse(Curve):
    r"""
    A class representing an Ellipse. \n

    If angle1 and angle2 are specified an EllipseArc between those angles is created.

    Attributes
    ----------
    center_x : float
        x-coordinate of the center point of th ellipse
    center_y : float
        y-coordinate of the center point of th ellipse
    center_z : float
        z-coordinate of the center point of th ellipse
    __radius_x : float
        radius in direction of the x-axis
    __radius_y : float
        radius in direction of the y-axis
    angle1 : float
        start angle of EllipseArc
    angle2 : float
        end angle of EllipseArc
    physical_tag : int, optional
        physiacal tag of the ellipse
    """

    __slots__ = ('center_x', 'center_y', 'center_z', '__radius_x',
                 '__radius_y', 'angle1', 'angle2')

    def __init__(self, center_x: float, center_y: float, center_z: float, radius_x: float, radius_y: float,
                 angle1: float = None, angle2: float = None, physical_tag: int = None) -> None:
        r"""
        The constructor of an Ellipse object. \n

        If angle1 and angle2 are specified an EllipseArc between those angles is created.

        Parameters
        ----------
        center_x : float
            x-coordinate of the center point of the ellipse
        center_y : float
            y-coordinate of the center point of the ellipse
        center_z : float
            z-coordinate of the center point of the ellipse
        radius_x : float
            radius in direction of the x-axis (has to be bigger than radius_y)
        radius_y : float
            radius in direction of the y-axis (has to be smaller than radius_x)
        angle1 : float, optional
            start angle of EllipseArc
        angle2 : float, optional
            end angle of EllipseArc

        physical_tag : int, optional
            physiacal tag of the ellipse, by default None

        Raises
        ------
        Exception
            Exception if radius == 0
        """
        self.center_x = center_x
        self.center_y = center_y
        self.center_z = center_z
        self.__radius_x = abs(radius_x)
        self.__radius_y = abs(radius_y)
        self.angle1 = angle1
        self.angle2 = angle2
        self.physical_tag = physical_tag

        # check radius
        if radius_x == 0:
            raise Exception("radius_x should be > 0.")
        if radius_y == 0:
            raise Exception("radius_y should be > 0.")
        if self.__radius_x < self.__radius_y:
            warnings.warn("Major radius radius_x should be larger than minor radius radius_y. "
                          "Rotate the shape or use a circle for radius_x = radius_y.")

        # create start and end point
        anglesGiven = False
        if self.angle1 is None or angle2 is None:
            angle1 = 0.
            if self.angle1 is not None or angle2 is not None:
                warnings.warn(
                    "Only one angle is set. Set two angles or none. The given angle currently is ignored.")
        else:
            anglesGiven = True
            angle1 = self.angle1
            angle2 = self.angle2
            x2, y2 = self.__calculate_point_on_ellipse(angle2)

        x1, y1 = self.__calculate_point_on_ellipse(angle1)

        start_point = Point(x1, y1, center_z)
        if anglesGiven:
            # otherwise end_point = start_point
            end_point = Point(x2, y2, center_z)
            super().__init__(start_point, end_point, physical_tag=physical_tag)
        else:
            super().__init__(start_point, physical_tag=physical_tag)

    @property
    def start_point(self) -> Point:
        """
        Returns start point of Ellips or EllipseArc. If Ellipse the start point is also the end point.

        Returns
        -------
        Point
            Start point of Ellips or EllipseArc

        Raises
        ------
        Warning
            If ellipse is added to gmsh start point might have other values.
        """
        if self.added is False:
            warnings.warn("The start point is not generated in gmsh yet. The values will be overwritten when added to "
                          "gmsh and might differ a bit.")
        return self.points[0]

    @property
    def end_point(self) -> Point:
        """
        Returns end point of Ellipse or EllipseArc. If ellipse the end point equals the start point.

        Returns
        -------
        Point
            End point of Ellips or EllipseArc

        Raises
        ------
        Warning
            If ellipse is added to gmsh end point might have other values.
        """
        if self.added is False:
            warnings.warn("The end point is not generated in gmsh yet. The values will be overwritten when added to "
                          "gmsh and might differ a bit.")
        return self.points[-1]

    @property
    def radius_x(self) -> float:
        """
        Getter property for radius_x

        Returns
        -------
        float
            radius in direction of the x-axis
        """
        return self.__radius_x

    @property
    def radius_y(self) -> float:
        """
        Getter property for radius_y

        Returns
        -------
        float
            radius in direction of the y-axis
        """
        return self.__radius_y

    @radius_x.setter
    def radius_x(self, radius_x) -> None:
        """
        Setter property for radius_x

        Parameters
        ----------
        radius_x : float
            radius in direction of the x-axis

        Raises
        ------
        Exception
            Exception if radius == 0
        """
        if radius_x == 0:
            raise Exception("radius_x should be > 0.")
        self.__radius_x = abs(radius_x)

    @radius_y.setter
    def radius_y(self, radius_y) -> None:
        """
        Setter property for radius_y

        Parameters
        ----------
        radius_y : float
            radius in direction of the y-axis

        Raises
        ------
        Exception
            Exception if radius == 0
        """
        if radius_y == 0:
            raise Exception("radius_y should be > 0.")
        self.__radius_y = abs(radius_y)

    def add_to_gmsh(self):
        """
        Adds an Ellipse or EllipseArc to gmsh.

        Raises
        ------
        Exception
            Exception if radius_x is smaller that radius_y
        Warning
            If Ellipse has already been added to gmsh
        """
        if self.added is False:
            if self.__radius_x < self.__radius_y:
                raise Exception(
                    "Major radius radius_x should be larger than minor radius radius_y")
            if self.angle1 is None or self.angle2 is None:
                self.tag = gmsh.model.occ.addEllipse(x=self.center_x, y=self.center_y, z=self.center_z,
                                                     r1=self.__radius_x, r2=self.__radius_y)
            else:
                self.tag = gmsh.model.occ.addEllipse(x=self.center_x, y=self.center_y, z=self.center_z,
                                                     r1=self.__radius_x, r2=self.__radius_y, angle1=self.angle1,
                                                     angle2=self.angle2)
            # add start and end point to self.points
            entities0D = gmsh.model.getEntities(0)
            for i, entity_0D in enumerate(entities0D):
                gmsh_point = gmsh.model.getValue(
                    entity_0D[0], entity_0D[1], parametricCoord=[])
                self.points[i].x = gmsh_point[0]
                self.points[i].y = gmsh_point[1]
                self.points[i].z = gmsh_point[2]
                self.points[i].tag = entity_0D[1]

        else:
            warnings.warn(
                "This Ellipse has already been added to gmsh. It will not be added again.")

    def __calculate_point_on_ellipse(self, angle: float) -> Point:
        """
        Returns the Point on the Ellipse at the given angle.

        Parameters
        ----------
        angle : float
            Angle in rad.

        Returns
        -------
        Point
            Point on Ellipse given an angle
        """
        # wrap angle from 0 to 2*pi
        angle = angle % (2 * math.pi)
        if angle < 0:
            angle = (angle + 2 * math.pi) % (2 * math.pi)
        # calculate theta
        theta = np.arctan(self.__radius_x * np.tan(angle) / self.__radius_y)
        if math.pi / 2 < angle < math.pi:
            theta += math.pi
        elif math.pi <= angle <= 3 * math.pi / 2:
            theta -= math.pi
        # calculate coordinates
        x = self.center_x + self.__radius_x * np.cos(theta)
        y = self.center_y + self.__radius_y * np.sin(theta)
        return x, y
