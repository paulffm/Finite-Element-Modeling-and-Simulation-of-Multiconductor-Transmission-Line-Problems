# coding=utf-8
"""Implementation of the class Yoke.

.. sectionauthor:: Bundschuh
"""
# pylint: disable=line-too-long

from typing import List, Union, Tuple
from dataclasses import dataclass
import numpy as np

from pyrit.geometry import Geo, Point, Line, CircleArc, Surface, NotAddedToGmshError

from pyrit import get_logger

logger = get_logger(__name__)


@dataclass
class Yoke(Geo):
    """Model of a Yoke with or without air gap in 2D. Meant for axialsymmetric geometries.

    Examples
    --------
    The easiest way to create this model is by

    >>> model = Yoke('yoke')

    This can then be used with the geometry module

    >>> from pyrit.geometry.Geometry import Geometry
    >>> geo = Geometry('yoke', show_gui=True)
    >>> with geo.geometry():
    >>>     model.add_to_gmsh()

    The parameters can be passed to the model with the kwargs:

    >>> model = Yoke('foil_winding', width=1.3, air_gap_pos=0.4)
    >>> model.corner_radius_inside = 0.05
    >>> model.corner_radius_outside_bot = 0.1
    """

    name: str = ''  #: Name of the foil winding
    origin_x: float = 0  #: x-coordinate of the origin
    origin_y: float = 0  #: y-coordinate of the origin
    origin_z: float = 0  #: z-coordinate of the origin
    anchor: Union[Geo.Anchor, Tuple[float]] = Geo.Anchor.WEST  #: Anchor of the foil winding.

    width: float = 1.  #: The width of the yoke (x direction).
    height: float = 2.  #: The height of the yoke (y direction).

    #: The thickness of the yoke. This is the default value that is used if the specific thicknesses for the sides
    #: (top, bot, inner, outer) are set to :py:const:`None`.
    thickness: float = 0.1
    thickness_bot: float = None  #: Thickness at the bottom side.
    thickness_top: float = None  #: Thickness at the top side.
    thickness_outer: float = None  #: Thickness at the outside, i.e. the right side.
    thickness_inner: float = None  #: Thickness in the middle, i.e. left side.

    #: The default value for the corner radius. It is used if `corner_radius_inside` or `corner_radius_outside` is set
    #: to :py:const:`None`. The corner radius can be 0. Then there is a sharp corner
    corner_radius: float = 0
    #: Corner radius for the inside corners. Used if a specific inner corner radius is set to :py:const:`None`.
    corner_radius_inside: float = None
    #: Corner radius for the outside corners. Used if a specific inner corner radius is set to :py:const:`None`.
    corner_radius_outside: float = None
    #: Corner radius for the inner (i.e. left) bottom corner at the inside.
    corner_radius_inside_bot_inner: float = None
    #: Corner radius for the outer (i.e. right) bottom corner at the inside.
    corner_radius_inside_bot_outer: float = None
    corner_radius_inside_top_inner: float = None  #: Corner radius for the inner (i.e. left) top corner at the inside.
    corner_radius_inside_top_outer: float = None  #: Corner radius for the outer (i.e. right) top corner at the inside.
    corner_radius_outside_bot: float = None  #: Corner radius for the bottom corner at the outside.
    corner_radius_outside_top: float = None  #: Corner radius for the top corner at the outside.

    air_gap_thickness: float = 0.05  #: The Thickness of the air gap.
    #: The position of the air gap. The middle of the air gap lies at `air_gap_pos` of the height. So with
    #: :code:`air_gap_pos=0.7` the air gap is at 70% of the height of the yoke.
    air_gap_pos: float = 0.5

    #: The mesh size for gmsh. This is the default mesh size taken for all points that are generated.
    mesh_size: float = None

    def __post_init__(self):
        super().__init__(self.name, self.origin_x, self.origin_y, self.origin_z, self.anchor)

        # protected variables
        self._inner_lines = None
        self._lines_axis = None
        self._outer_lines = None
        self._surface = None

    def scale(self, scale: float):
        parameters = ('width', 'height', 'thickness', 'thickness_bot', 'thickness_outer', 'thickness_top',
                      'thickness_inner', 'corner_radius', 'corner_radius_inside', 'corner_radius_outside',
                      'corner_radius_inside_bot_inner', 'corner_radius_inside_bot_outer',
                      'corner_radius_inside_top_inner', 'corner_radius_inside_top_outer', 'corner_radius_outside_bot',
                      'corner_radius_outside_top', 'air_gap_thickness')
        for parameter in parameters:
            if self.__getattribute__(parameter) is not None:
                self.__setattr__(parameter, scale * self.__getattribute__(parameter))

    def scale_x(self, scale: float):
        """Scale the model in x-direction.

        Parameters
        ----------
        scale : float
            The factor to scale.
        """
        if self.width is not None:
            self.width = scale * self.width
        for parameter in ('thickness_inner', 'thickness_outer'):
            self.__setattr__(parameter, scale * self.get_thickness(parameter))

    def scale_y(self, scale: float):
        """Scale the model in y-direction.

        Parameters
        ----------
        scale : float
            The factor to scale.
        """
        if self.height is not None:
            self.height = scale * self.height
        if self.air_gap_thickness is not None:
            self.air_gap_thickness = scale * self.air_gap_thickness
        for parameter in ('thickness_bot', 'thickness_top'):
            self.__setattr__(parameter, scale * self.get_thickness(parameter))

    def get_corner_radius(self, key: str) -> float:
        """Determines the corner radius of a specific corner.

        Parameters
        ----------
        key : str
            The identifier of a parameter. See :py:attr:`possible_parameters`.

        Returns
        -------
        corner_radius : float
            The corner radius of the corner
        """
        possible_parameters = ('corner_radius', 'corner_radius_inside', 'corner_radius_outside',
                               'corner_radius_inside_bot_inner', 'corner_radius_inside_bot_outer',
                               'corner_radius_inside_top_inner', 'corner_radius_inside_top_outer',
                               'corner_radius_outside_bot', 'corner_radius_outside_top')
        if key not in possible_parameters:
            raise KeyError(f"The key '{key}' is not in the possible parameters.")
        if self.__getattribute__(key) is None:
            if 'inside' in key:
                if self.corner_radius_inside is None:
                    return self.corner_radius
                return self.corner_radius_inside
            if 'outside' in key:
                if self.corner_radius_outside is None:
                    return self.corner_radius
                return self.corner_radius_outside
            return None
        return self.__getattribute__(key)

    def get_thickness(self, key: str) -> float:
        """Returns th thickness of a specific side.

        Parameters
        ----------
        key : str
            Identifier of the side. See :py:attr:`possible_parameters`.

        Returns
        -------
        thickness : float
            The thickness of this side.
        """
        possible_parameters = ('thickness', 'thickness_bot', 'thickness_outer', 'thickness_top', 'thickness_inner')
        if key not in possible_parameters:
            raise KeyError(f"The key '{key}' is not in the possible parameters.")
        if self.__getattribute__(key) is None:
            return self.thickness
        return self.__getattribute__(key)

    @staticmethod
    def _help_draw(p: Point, radius: float, fun, dx: float = 0, dy: float = 0):
        """Helps to generate the geometry."""
        if (dx == 0 and dy == 0) or (dx != 0 and dy != 0):
            raise Exception

        lines = []
        if radius == 0:
            last_point = fun(p.x + dx, p.y + dy)
            lines.append(Line(p, last_point))
            return lines, last_point, None

        ddx, ddy = 0, 0
        if dx != 0:
            ddy = dx
            ddx = 0
        if dy != 0:
            ddx = -dy
            ddy = 0
        start_point = fun(p.x + dx - np.sign(dx) * radius, p.y + dy - np.sign(dy) * radius)
        end_point = fun(p.x + dx + np.sign(ddx) * radius, p.y + dy + np.sign(ddy) * radius)
        center_point = fun(p.x + dx - np.sign(dx) * radius + np.sign(ddx) * radius,
                           p.y + dy - np.sign(dy) * radius + np.sign(ddy) * radius)
        arc = CircleArc(start_point, end_point, center_point)
        lines.append(Line(p, arc.start_point))
        lines.append(arc)
        return lines, arc.end_point, arc.center_point

    def determine_coordinates(self, anchor: Union[Geo.Anchor, Tuple[float]]) -> np.ndarray:
        """Computes the coordinates of an anchor.

        Parameters
        ----------
        anchor : Union[Anchor, Tuple[float]]
            An anchor or tuple of the relative coordinates. (0,0) is the south-west corner (1,1) is the north-east
            corner.

        Returns
        -------
        coordinates : np.ndarray
            The coordinates of the anchor.
        """
        coordinates = self.get_origin()
        if isinstance(self.anchor, self.Anchor):
            self_anchor = np.array(self.anchor.value)
        else:
            self_anchor = np.array(self.anchor)
        anchor = np.array(anchor.value) - self_anchor
        coordinates[:2] = coordinates[:2] + np.multiply(anchor, [self.width, self.height])
        return coordinates

    # noinspection PyAttributeOutsideInit
    def add_to_gmsh(self, **kwargs):
        height = self.height
        width = self.width

        r_in_bot_in = self.get_corner_radius('corner_radius_inside_bot_inner')
        r_in_bot_out = self.get_corner_radius('corner_radius_inside_bot_outer')
        r_in_top_in = self.get_corner_radius('corner_radius_inside_top_inner')
        r_in_top_out = self.get_corner_radius('corner_radius_inside_top_outer')
        r_out_bot = self.get_corner_radius('corner_radius_outside_bot')
        r_out_top = self.get_corner_radius('corner_radius_outside_top')

        agt = self.air_gap_thickness
        agp = self.air_gap_pos

        gap_y_bot = agp * height - 0.5 * agt
        gap_y_top = agp * height + 0.5 * agt

        t_bot = self.get_thickness('thickness_bot')
        t_out = self.get_thickness('thickness_outer')
        t_top = self.get_thickness('thickness_top')
        t_in = self.get_thickness('thickness_inner')

        mesh_size = self.mesh_size

        south_west = self.determine_coordinates(self.Anchor.SOUTH_WEST)

        z = south_west[2]

        def add_point_offset(x, y):
            return Point(x + south_west[0], y + south_west[1], z, mesh_size=mesh_size)

        def add_point(x, y):
            return Point(x, y, z, mesh_size=mesh_size)

        point_gap_in_top = point_gap_out_top = point_gap_in_bot = point_gap_out_bot = None

        # region Lines outside
        lines_outside = []
        del_points = []

        first_point = add_point_offset(0, 0)
        last_point = first_point

        lines, last_point, del_point = self._help_draw(last_point, r_out_bot, add_point, dx=width)
        lines_outside.extend(lines)
        del_points.append(del_point)

        if r_out_bot == 0:
            dy = height
        else:
            dy = height - r_out_bot
        lines, last_point, del_point = self._help_draw(last_point, r_out_top, add_point, dy=dy)
        lines_outside.extend(lines)
        del_points.append(del_point)

        # lines_outside.append(Line(last_point, add_point(0, height)))

        if agt == 0:  # No air gap
            points_tmp = [last_point, add_point_offset(0, height), first_point]
            tmp = -1
        else:
            point_gap_out_top = add_point_offset(0, gap_y_top)
            point_gap_out_bot = add_point_offset(0, gap_y_bot)
            points_tmp = [last_point, add_point_offset(0, height), point_gap_out_top, point_gap_out_bot, first_point]
            tmp = -3
        for k in range(len(points_tmp) - 1):
            lines_outside.append(Line(points_tmp[k], points_tmp[k + 1]))

        self._outer_lines = lines_outside[:tmp]
        self._lines_axis = lines_outside[tmp:]
        # endregion

        # region Lines inside
        first_point_inside = None
        lines_inside = []

        if agt == 0:  # No air gap
            if r_in_bot_in == 0:
                first_point_inside = add_point_offset(t_in, t_bot)
                last_point = first_point_inside
            else:
                arc = CircleArc(add_point_offset(t_in, t_bot + r_in_bot_in),
                                add_point_offset(t_in + r_in_bot_in, t_bot),
                                add_point_offset(t_in + r_in_bot_in, t_bot + r_in_bot_in))
                lines_inside.append(arc)
                first_point_inside = arc.start_point
                last_point = arc.end_point
                del_points.append(arc.center_point)
        else:  # Air gap
            point_gap_in_bot = add_point_offset(t_in, gap_y_bot)
            lines, last_point, del_point = self._help_draw(point_gap_in_bot, r_in_bot_in, add_point,
                                                           dy=-gap_y_bot + t_bot)
            lines_inside.extend(lines)
            del_points.append(del_point)

        # Bottom line and bottom right arc
        if r_in_bot_in == 0:
            dx = width - t_in - t_out
        else:
            dx = width - t_in - t_out - r_in_bot_in
        lines, last_point, del_point = self._help_draw(last_point, r_in_bot_out, add_point, dx=dx)
        lines_inside.extend(lines)
        del_points.append(del_point)

        # Right line and top right arc
        if r_in_bot_out == 0:
            dy = height - t_bot - t_top
        else:
            dy = height - t_bot - t_top - r_in_bot_out
        lines, last_point, del_point = self._help_draw(last_point, r_in_top_out, add_point, dy=dy)
        lines_inside.extend(lines)
        del_points.append(del_point)

        # Top line and top left arc
        if r_in_top_out == 0:
            dx = - width + t_in + t_out
        else:
            dx = - width + t_in + t_out + r_in_top_out
        lines, last_point, del_point = self._help_draw(last_point, r_in_top_in, add_point, dx=dx)
        lines_inside.extend(lines)
        del_points.append(del_point)

        if agt == 0:  # No air gap
            lines_inside.append(Line(last_point, first_point_inside))
        else:
            point_gap_in_top = add_point_offset(t_in, gap_y_top)
            lines_inside.append(Line(last_point, point_gap_in_top))
        # endregion

        for line in lines_outside:
            line.add_to_gmsh()
        for line in lines_inside:
            line.add_to_gmsh()

        if agt == 0:  # No air gap
            surf = Surface(lines_outside, lines_inside)
            surf.add_to_gmsh()
            self._inner_lines = lines_inside
        else:
            line_gap_top = Line(point_gap_in_top, point_gap_out_top)
            line_gap_bot = Line(point_gap_in_bot, point_gap_out_bot)

            surf = Surface(lines_outside[:-2] + lines_outside[-1:] + [line_gap_top, line_gap_bot] + lines_inside)
            surf.add_to_gmsh()
            self._inner_lines = lines_inside + [line_gap_top, line_gap_bot, lines_outside[-2]]
        for p in del_points:
            if p is not None:
                p.remove()

        self._surface = surf

    @property
    def inner_lines(self) -> List[Line]:
        """List of the lines at the inside of the yoke."""
        if self._inner_lines is None:
            raise NotAddedToGmshError('')
        return self._inner_lines

    @property
    def outer_lines(self) -> List[Line]:
        """List of lines at the outside of the yoke. Lines on the axis (i.e. left side) are not in this list."""
        if self._outer_lines is None:
            raise NotAddedToGmshError('')
        return self._outer_lines

    @property
    def axis_lines(self):
        """List of lines at the axis (i.e. left side)."""
        if self._lines_axis is None:
            raise NotAddedToGmshError('')
        return self._lines_axis

    @property
    def surface(self) -> Surface:
        """The surface of the yoke."""
        if self._surface is None:
            raise NotAddedToGmshError('')
        return self._surface
