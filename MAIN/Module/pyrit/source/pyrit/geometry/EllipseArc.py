# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 10:49:49 2021

.. sectionauthor:: Koray
"""

import gmsh


from . import Curve, Point


class EllipseArc(Curve):
    """Class representing an ellipse"""

    def __init__(self, start_point: Point, end_point: Point, center_point: Point,
                 major_point: Point, physical_tag: int = None):
        """
        Constructor of the class EllipseArc.

        Parameters
        ----------
        start_point : Point
            Starting point of the ellipse arc.
        end_point : Point
            Ending point of the ellipse arc.
        center_point : Point
            Center point of the ellipse arc.
        major_point : Point
            Major point of the ellipse arc.
        physical_tag : int, optional
            Physical tag of the entity. The default is None.

        Returns
        -------
        None.

        """
        super().__init__(start_point, end_point, center_point, major_point, physical_tag=physical_tag)

    @property
    def start_point(self) -> Point:
        """
        Get the start point of the ellipse arc.

        Returns
        -------
        Point
            start point.

        """
        return self.points[0]

    @property
    def end_point(self) -> Point:
        """
        Get the end point of the ellipse arc.

        Returns
        -------
        Point
            end point.

        """
        return self.points[1]

    @property
    def center_point(self) -> Point:
        """
        Get the center point of the ellipse arc.

        Returns
        -------
        Point
            center point.

        """
        return self.points[2]

    @property
    def major_point(self) -> Point:
        """
        Get the major point of the ellipse arc.

        Returns
        -------
        Point
            Major point.

        """
        return self.points[3]

    def add_to_gmsh(self):
        """
        Adds Ellipse Arc to GMSH.

        Returns
        -------
        None.

        """
        for pt in self.points:
            if pt.added is False:
                pt.add_to_gmsh()
        self.tag = gmsh.model.occ.addEllipseArc(
            startTag=self.start_point.tag, centerTag=self.center_point.tag,
            majorTag=self.major_point.tag, endTag=self.end_point.tag)
