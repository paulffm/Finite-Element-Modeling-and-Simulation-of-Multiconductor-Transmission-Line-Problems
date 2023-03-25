# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 2021

Contains the class CircleArc

.. sectionauthor:: menzenbach
"""

import gmsh


from . import Point, Curve


class CircleArc(Curve):
    """A class representing a CircleArc."""

    def __init__(self, start_point: Point, end_point: Point, center_point: Point, physical_tag: int = None) -> None:
        """
        Constructor of the class CircleArc

        Parameters
        ----------
        start_point : Point
            Start point of the Arc
        end_point : Point
            End point of the Arc
        center_point : Point
            Center Point of the Arc
        physical_tag : int, optional
            Physical tag of the entity, by default None
        """
        super().__init__(start_point, end_point, center_point, physical_tag=physical_tag)

    @property
    def start_point(self) -> Point:
        """
        Get the start point of the circle arc.

        Returns
        -------
        Point
            Start point.

        """
        return self.points[0]

    @property
    def end_point(self) -> Point:
        """
        Get the end point of the circle arc.

        Returns
        -------
        Point
            End point.

        """
        return self.points[1]

    @property
    def center_point(self) -> Point:
        """
        Get the center point of the circle arc.

        Returns
        -------
        Point
            Center point.

        """
        return self.points[2]

    def add_to_gmsh(self):
        for pt in self.points:
            if pt.added is False:
                pt.add_to_gmsh()
        self.tag = gmsh.model.occ.addCircleArc(
            startTag=self.start_point.tag, centerTag=self.center_point.tag, endTag=self.end_point.tag)
