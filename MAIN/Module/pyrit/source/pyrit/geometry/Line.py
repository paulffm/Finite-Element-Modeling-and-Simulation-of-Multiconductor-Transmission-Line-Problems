# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 2021

Contains the class Line.

.. sectionauthor:: leissner
"""

import gmsh

from . import Curve, Point


class Line(Curve):
    """
    A class representing a Line.

    Methods
    -------
    add_in_GMSH()
        add the Line in GMSH
    """

    def __init__(self, point1: Point, point2: Point, physical_tag: int = None):
        """
        Initialize a Line.

        Parameters
        ----------
        point1 : Point
            First point of the line.
        point2 : Point
            Second point of the line.
        physical_tag : int, optional
            Physical tag of the entity. The default is None.

        Returns
        -------
        None.

        """
        super().__init__(point1, point2, physical_tag=physical_tag)

    def __repr__(self):
        return f"Line with tag '{self._tag}' between Points '{self.points[0]}' and '{self.points[1]}'"

    def add_to_gmsh(self):
        """
        Add Line in GMSH.

        Returns
        -------
        None.

        """
        for pt in self.points:
            if pt.added is False:
                pt.add_to_gmsh()
        self.tag = gmsh.model.occ.addLine(self.points[0].tag, self.points[1].tag)
