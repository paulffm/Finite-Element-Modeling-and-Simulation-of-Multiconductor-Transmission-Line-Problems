# coding=utf-8
"""
Class FoilWinding

.. sectionauthor:: Bundschuh
"""

from typing import NoReturn, List, Union, Tuple
from dataclasses import dataclass
import numpy as np

from pyrit.geometry import Geo, Point, Line, Surface, NotAddedToGmshError

from pyrit import get_logger

logger = get_logger(__name__)


@dataclass
class FoilWinding(Geo):
    """Model of foil windings in a 2D geometry.

    There are optional paramters. Most of them are self-explanatory. The others are explained here:

        - anchor
            This is the position of the foil winding that lies on the origin coordinates. It can be a predefined
            Anchor (:py:class:`Anchor`) or a tuple of two numbers. (0,0) is the south-west corner (1,1) is
            the north-east corner.
        - foil_insulation_distribution
            Determines how the cinductor and insulator of a foil are arranged. This is a number in :math:`(0,1)`, where
            0 means that the conductor is as much left as possible and 1 means that it is as much richt as possible. The
            value of 0.5 for exapmle means that on each side of the conductor, there is an insulation layer of the same
            thickness.

    Note that the layout of the windings can be determined **either** by :py:attr:`foil_thickness` and
    :py:attr:`insulation_thickness` (default) **or** :py:attr:`fill_factor` and :py:attr:`width`. After creation, the
    other values will be set appropriately.


    Examples
    --------
    The easiest way to create this model is by

    >>> model = FoilWinding('foil_winding')

    This can then be used with the geometry module

    >>> from pyrit.geometry.Geometry import Geometry
    >>> geo = Geometry('foil_winding', show_gui=True)
    >>> with geo.geometry():
    >>>     model.add_to_gmsh()
    """

    name: str = ''  #: Name of the foil winding
    origin_x: float = 0  #: x-coordinate of the origin
    origin_y: float = 0  #: y-coordinate of the origin
    origin_z: float = 0  #: z-coordinate of the origin
    anchor: Union[Geo.Anchor, Tuple[float]] = Geo.Anchor.SOUTH_WEST  #: Anchor of the foil winding.
    height: float = 1.  #: Height of the foil winding
    turns: int = 10  #: Number of turns of the foil winding
    angle: float = 0  #: Angle of the foil winding
    foil_insulation_distribution: float = 0.5  #: How foil and insulation are arranged. See class documentation.
    mesh_size: float = 0.1  #: Mesh size inside the foil_winding region

    foil_thickness: float = 0.1  #: Thickness of the foil
    insulation_thickness: float = 0.01  #: Thickness of the insulation

    fill_factor: float = None  #: The fill factor
    width: float = None  #: The total width

    def __post_init__(self):
        super().__init__(self.name, self.origin_x, self.origin_y, self.origin_z, self.anchor)
        self.use_thicknesses = True
        if self.fill_factor is not None and self.width is not None:  # both are set
            logger.info("Using the fill factor and width.")
            self.foil_thickness = self.width / self.turns * self.fill_factor
            self.insulation_thickness = self.width / self.turns * (1 - self.fill_factor)
            self.use_thicknesses = False
        else:
            if self.fill_factor is None and self.width is None:  # None is set
                logger.info("Using the thickness of foil and insulation.")
            else:
                logger.warning(
                    "Only one of fill_factor and width were given but both are needed. Using the default settings.")
            self.width = self.turns * (self.insulation_thickness + self.foil_thickness)
            self.fill_factor = self.foil_thickness / (self.foil_thickness + self.insulation_thickness)

        # Check foil_insulation_distribution
        if not 0 < self.foil_insulation_distribution < 1:
            logger.warning('The foil_insulation_distribution should be between 0 and 1 (exclusive). '
                           'This may result in strange geometries.')

        # protected variables
        self._points = None
        self._lines = None
        self._surfaces = None

    def determine_coordinates(self, anchor: Union[Geo.Anchor, Tuple[float]]) -> np.ndarray:
        coordinates = self.get_origin()
        if isinstance(self.anchor, self.Anchor):
            self_anchor = np.array(self.anchor.value)
        else:
            self_anchor = np.array(self.anchor)
        anchor = np.array(anchor.value) - self_anchor
        unit_width = self.width * np.array([np.cos(self.angle), np.sin(self.angle)])
        unit_height = self.height * np.array([-1 * np.sin(self.angle), np.cos(self.angle)])
        coordinates[:2] = coordinates[:2] + anchor[0] * unit_width + anchor[1] * unit_height
        return coordinates

    def scale(self, scale: float) -> NoReturn:
        self.scale_x(scale)
        self.scale_y(scale)

    def scale_x(self, scale: float) -> NoReturn:
        """Scale the model in x direction.

        Parameters
        ----------
        scale : float
            The factor to scale in x direction.
        """
        self.width = scale * self.width
        self.foil_thickness = scale * self.foil_thickness
        self.insulation_thickness = scale * self.insulation_thickness

    def scale_y(self, scale: float) -> NoReturn:
        """Scale the model in y direction.

        Parameters
        ----------
        scale : float
            The factor to scale in y direction.
        """
        self.height = scale * self.height

    def add_to_gmsh(self, **kwargs):
        mesh_size = self.mesh_size

        south_west = self.determine_coordinates(self.Anchor.SOUTH_WEST)

        wf = self.foil_thickness
        wi = self.insulation_thickness

        coords_bot = south_west[0:2]
        coords_top = self.determine_coordinates(self.Anchor.NORTH_WEST)[0:2]
        z = south_west[2]

        def shift_width(coords, shift):
            return coords + shift * np.array([np.cos(self.angle), np.sin(self.angle)])

        points_bot = [Point(*coords_bot, z, mesh_size=mesh_size)]
        points_top = [Point(*coords_top, z, mesh_size=mesh_size)]
        coords_bot = shift_width(coords_bot, self.foil_insulation_distribution * wi)
        coords_top = shift_width(coords_top, self.foil_insulation_distribution * wi)
        for _ in range(self.turns):
            points_bot.append(Point(*coords_bot, z, mesh_size=mesh_size))
            coords_bot = shift_width(coords_bot, wf)
            points_bot.append(Point(*coords_bot, z, mesh_size=mesh_size))

            points_top.append(Point(*coords_top, z, mesh_size=mesh_size))
            coords_top = shift_width(coords_top, wf)
            points_top.append(Point(*coords_top, z, mesh_size=mesh_size))

            coords_bot = shift_width(coords_bot, wi)
            coords_top = shift_width(coords_top, wi)

        coords_bot = shift_width(coords_bot, - (1 - self.foil_insulation_distribution) * wi)
        coords_top = shift_width(coords_top, - (1 - self.foil_insulation_distribution) * wi)

        points_bot.append(Point(*coords_bot, z, mesh_size=mesh_size))
        points_top.append(Point(*coords_top, z, mesh_size=mesh_size))

        lines_top = []
        lines_bot = []
        lines_vertical = [Line(points_top[0], points_bot[0])]
        lines_vertical[0].add_to_gmsh()
        for k in range(len(points_top) - 1):
            lines_top.append(Line(points_top[k], points_top[k + 1]))
            lines_top[k].add_to_gmsh()
            lines_bot.append(Line(points_bot[k], points_bot[k + 1]))
            lines_bot[k].add_to_gmsh()
            lines_vertical.append(Line(points_top[k + 1], points_bot[k + 1]))
            lines_vertical[-1].add_to_gmsh()

        foil_surfaces = []
        for k in range(self.turns):
            lines = [lines_top[1 + 2 * k], lines_vertical[2 + 2 * k], lines_bot[1 + 2 * k], lines_vertical[1 + 2 * k]]
            surf = Surface(lines)
            surf.add_to_gmsh()
            foil_surfaces.append(surf)

        iso_surfaces = []
        for k in range(self.turns + 1):
            lines = [lines_top[2 * k], lines_vertical[1 + 2 * k], lines_bot[2 * k], lines_vertical[2 * k]]
            surf = Surface(lines)
            surf.add_to_gmsh()
            iso_surfaces.append(surf)

        # noinspection PyAttributeOutsideInit
        self._points = [points_bot, points_top]
        # noinspection PyAttributeOutsideInit
        self._lines = [lines_bot, lines_top, lines_vertical]
        # noinspection PyAttributeOutsideInit
        self._surfaces = [foil_surfaces, iso_surfaces]

    @property
    def points(self) -> List[List[Point]]:
        """Returns all Points of the model.

        Returns
        -------
        out : List[List[Point]]
            The first entry is a list of the bottom points and the second entry a list ot the top points.
        """
        if self._points is None:
            raise NotAddedToGmshError('')
        return self._points

    @property
    def lines(self) -> List[List[Line]]:
        """Returns all lines of the model.

        Returns
        -------
        out : List[List[Line]]
            The first entry is a list of the bottom lines, the second entry of the top lines and the third entry of the
            vertical lines.
        """
        if self._lines is None:
            raise NotAddedToGmshError('')
        return self._lines

    @property
    def foil_surfaces(self) -> List[Surface]:
        """A list of all foil surfaces.

        Returns
        -------
        out : List[Surface]
            All surfaces that represent a foil.
        """
        if self._surfaces is None:
            raise NotAddedToGmshError('')
        return self._surfaces[0]

    @property
    def insulation_surfaces(self) -> List[Surface]:
        """A list of all insulation surfaces.

        Returns
        -------
        out : List[Surface]
            All surfaces that represent an insulation
        """
        if self._surfaces is None:
            raise NotAddedToGmshError('')
        return self._surfaces[1]

    @property
    def outer_lines(self) -> List[Line]:
        """A list of all outer lines.

        This list can be used for the definition of another surface around the foil winding.

        Returns
        -------
        out : List[Line]
            List of all outer lines.
        """
        if self._lines is None:
            raise NotAddedToGmshError('')
        return self._lines[0] + [self._lines[2][-1]] + [self._lines[1][-k] for k in
                                                        range(1, len(self._lines[1]) + 1)] + [self._lines[2][0]]
