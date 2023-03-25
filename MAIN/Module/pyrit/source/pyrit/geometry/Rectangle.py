# coding=utf-8
"""
A rectangle

.. sectionauthor:: bundschuh
"""

from typing import Type
import gmsh
import numpy as np

from . import Surface, Point, Curve


class Rectangle(Surface):
    """Class representing a surface."""

    __slots__ = ('origin_x', 'origin_y', 'origin_z', '__dx', '__dy', '__corner_radius')

    def __init__(self, origin_x: float, origin_y: float, origin_z: float, dx: float, dy: float,
                 corner_radius: float = None, physical_tag: int = None):
        super().__init__([], physical_tag=physical_tag)
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.origin_z = origin_z

        self.dx = dx
        self.dy = dy
        self.corner_radius = corner_radius

    # region Properties

    @property
    def dx(self) -> float:
        """Returns the x width of the rectangle."""
        return self.__dx

    @dx.setter
    def dx(self, dx: float) -> None:
        """Set the x width of the rectangle."""
        if dx > 0:
            self.__dx = dx
        else:
            raise ValueError(f"dx has to be a positive number, but was {dx}")

    @property
    def dy(self) -> float:
        """Returns the y width of the rectangle."""
        return self.__dy

    @dy.setter
    def dy(self, dy: float) -> None:
        """Set the y width of the rectangle."""
        if dy > 0:
            self.__dy = dy
        else:
            raise ValueError(f"dy has to be a positive number, but was {dy}")

    @property
    def corner_radius(self) -> float:
        """Get the corner radius of the rectangle."""
        return self.__corner_radius

    @corner_radius.setter
    def corner_radius(self, corner_radius: float) -> None:
        """Set the corner radius of the rectangle."""
        if corner_radius is None:
            self.__corner_radius = corner_radius
        elif corner_radius > 0:
            if corner_radius < self.dx / 2 and corner_radius < self.dy / 2:
                self.__corner_radius = corner_radius
            else:
                raise ValueError("corner_radius has to be smaller than dx/2 and dy/2.")
        else:
            raise ValueError(f"corner_radius has to be a positive number, but was {corner_radius}")

    @property
    def has_rounded_corners(self) -> bool:
        """Returns True of the rectangle has rounded corners."""
        return self.corner_radius is not None

    @property
    def curve_south(self) -> Type[Curve]:
        """Return the south curve."""
        return self.ext_curves[0]

    @property
    def curve_east(self) -> Type[Curve]:
        """Return the east curve."""
        if self.has_rounded_corners:
            return self.ext_curves[2]
        return self.ext_curves[1]

    @property
    def curve_north(self) -> Type[Curve]:
        """Return the north curve."""
        if self.has_rounded_corners:
            return self.ext_curves[4]
        return self.ext_curves[2]

    @property
    def curve_west(self) -> Type[Curve]:
        """Return the west curve."""
        if self.has_rounded_corners:
            return self.ext_curves[6]
        return self.ext_curves[3]

    @property
    def curve_south_west(self) -> Type[Curve]:
        """Return the south west curve."""
        if self.has_rounded_corners:
            return self.ext_curves[7]
        return None

    @property
    def curve_south_east(self) -> Type[Curve]:
        """Return the south east curve."""
        if self.has_rounded_corners:
            return self.ext_curves[1]
        return None

    @property
    def curve_north_east(self) -> Type[Curve]:
        """Return the north east curve."""
        if self.has_rounded_corners:
            return self.ext_curves[3]
        return None

    @property
    def curve_north_west(self) -> Type[Curve]:
        """Return the north west curve."""
        if self.has_rounded_corners:
            return self.ext_curves[5]
        return None

    @property
    def point_south_west(self) -> Point:
        """Return the south west point."""
        return self.curve_south.points[0]

    @property
    def point_south_east(self) -> Point:
        """Return the south east point."""
        return self.curve_south.points[1]

    @property
    def point_north_east(self) -> Point:
        """Return the north east point."""
        return self.curve_north.points[0]

    @property
    def point_north_west(self) -> Point:
        """Return the north west point."""
        return self.curve_north.points[1]

    @property
    def point_west_south(self) -> Point:
        """Return the west south point."""
        return self.curve_west.points[1]

    @property
    def point_west_north(self) -> Point:
        """Return the west north point."""
        return self.curve_west.points[0]

    @property
    def point_east_south(self) -> Point:
        """Return the east south point."""
        return self.curve_east.points[0]

    @property
    def point_east_north(self) -> Point:
        """Return the east north point."""
        return self.curve_east.points[1]

    # endregion

    def add_to_gmsh(self):
        tags_points_before = np.array([point[1] for point in gmsh.model.occ.getEntities(0)])
        tags_curves_before = np.array([curve[1] for curve in gmsh.model.occ.getEntities(1)])
        if self.has_rounded_corners:
            self.tag = gmsh.model.occ.addRectangle(self.origin_x, self.origin_y, self.origin_z, self.__dx, self.__dy,
                                                   roundedRadius=self.corner_radius)
        else:
            self.tag = gmsh.model.occ.addRectangle(self.origin_x, self.origin_y, self.origin_z, self.__dx, self.__dy)
        tags_points_after = np.array([point[1] for point in gmsh.model.occ.getEntities(0)])
        tags_curves_after = np.array([curve[1] for curve in gmsh.model.occ.getEntities(1)])

        points_diff = np.setdiff1d(tags_points_after, tags_points_before)
        curves_diff = np.setdiff1d(tags_curves_after, tags_curves_before)

        if self.has_rounded_corners:
            x_values = np.array([self.__corner_radius, self.__dx - self.__corner_radius, self.__dx, self.__dx,
                                 self.__dx - self.__corner_radius, self.__corner_radius, 0, 0]) + self.origin_x
            y_values = np.array([0, 0, self.corner_radius, self.__dy - self.__corner_radius, self.__dy, self.__dy,
                                 self.__dy - self.__corner_radius, self.__corner_radius]) + self.origin_y
        else:
            x_values = np.array([0, self.__dx, self.__dx, 0]) + self.origin_x
            y_values = np.array([0, 0, self.__dy, self.__dy]) + self.origin_y

        points = [Point(x, y, self.origin_z) for x, y in zip(x_values, y_values)]
        for k, point in enumerate(points):
            point.tag = points_diff[k]

        curves = [Curve(points[k], points[k + 1]) for k, _ in enumerate(points[:-1])]
        curves.append(Curve(points[-1], points[0]))
        for k, curve in enumerate(curves):
            curve.tag = curves_diff[k]

        self.ext_curves = curves
        return self.tag
