# coding=utf-8
"""Foil windings with a surrounding.

.. sectionauthor:: Bundschuh
"""

from typing import NoReturn, List, Union, Tuple
from dataclasses import dataclass
import numpy as np
import gmsh

from pyrit import get_logger

from pyrit.geometry import Geo, Point, Line, Surface, CircleArc, Curve
from .FoilWinding import FoilWinding

logger = get_logger(__name__)


# noinspection PyAttributeOutsideInit
@dataclass
class FoilWindingWithSurrounding(Geo):
    """Model of foil windings in a 2D geometry with a surrounding.

    This class gets an instance of :py:class:`FoilWinding` and settings on how the surrounding looks like. The two
    options are a rectangular bounding box (default) or a round bounding box. If :py:attr:`radius_factor` is given, the
    round bounding box is used.

    To hold the size of the bounding box flexible, it conforms with the size of the foils windings. For that reason,
    this class works with factors. For the rectangular box there are factors for top, bottom and the right side. They
    indicate how much larger than the height and the outer radius, respectively, the bounding box is. For the round
    bounding box, there is only a factor for the radius.
    """

    foil_winding: FoilWinding  #: An instance of FoilWinding
    name: str = ''  #: Name of the model
    factor_bottom: float = 3  #: The factor for the bottom side
    factor_top: float = 3  #: The factor for the top side
    factor_right: float = 3  #: The factor for the right side
    factor: float = None  #: Joint factor for rectangular bounding box.

    radius_factor: float = None  #: The factor for the radius

    def __post_init__(self):
        origin_y = self.foil_winding.determine_coordinates(Geo.Anchor.WEST)[1]
        super().__init__(self.name, 0, origin_y, 0, self.Anchor.WEST)
        if self.factor is not None and self.radius_factor is not None:
            logger.warning("factor and radius_factor given. Using tha latter")
            self.factor = None

        if self.factor is not None:
            logger.info("Using a rectangular bounding box.")
            self.factor_bottom = self.factor
            self.factor_right = self.factor
            self.factor_top = self.factor

        if self.radius_factor is not None:
            logger.info("Using a round bounding box.")
            if self.radius_factor <= 0:
                raise ValueError(f"radius factor needs to be greater than 0, but is {self.radius_factor=}")

        self.tag_distance = None
        self.tag_threshold = None

        self._curves = None
        self._outer_curves = None
        self._surrounding = None

    def scale(self, scale: float) -> NoReturn:
        self.foil_winding.scale(scale)

    def scale_x(self, scale: float) -> NoReturn:
        """Scale the model in x direction.

        Parameters
        ----------
        scale : float
            The factor to scale in x direction.
        """
        self.foil_winding.scale_x(scale)

    def scale_y(self, scale: float) -> NoReturn:
        """Scale the model in y direction.

        Parameters
        ----------
        scale : float
            The factor to scale in y direction.
        """
        self.foil_winding.scale_y(scale)

    @property
    def height(self):
        """Height of the model."""
        if self.radius_factor is None:
            return (1 + self.factor_bottom + self.factor_top) * self.foil_winding.height
        return 2 * self.radius

    @property
    def width(self):
        """Width of the model."""
        if self.radius_factor is None:
            return self.foil_winding.east[0] * (1 + self.factor_right)
        return self.radius

    @property
    def radius(self):
        """Radius of the model. None if the rectangular bounding box is used."""
        if self.radius_factor is None:
            return None
        base_radius = np.sqrt(self.foil_winding.east[0] ** 2 + (0.5 * self.foil_winding.height ** 2))
        return base_radius + base_radius * self.radius_factor

    def determine_coordinates(self, anchor: Union[Geo.Anchor, Tuple[float]]) -> np.ndarray:
        coordinates = self.get_origin()
        if isinstance(self.anchor, self.Anchor):
            self_anchor = np.array(self.anchor.value)
        else:
            self_anchor = np.array(self.anchor)
        anchor = np.array(anchor.value) - self_anchor
        coordinates[:2] = coordinates[:2] + np.multiply(anchor, [self.width, self.height])
        return coordinates

    def add_to_gmsh(self, **kwargs):
        """Adds the Entities of the `Geo` to gmsh.

        When this method is called, gmsh needs to be initialized. So at least

        .. code-block:: python

            gmsh.initialize(sys.argv)

        has to be called beforehand. However, it is recommended to use the context manager
        :py:class:`~pyrit.geometry.Geometry.Geometry` from the class `Geometry`.

        Parameters
        ----------
        kwargs :
            Additional keyword arguments. Options are:

                mesh_size_inner : float
                    The mesh size in the middle of the coil (at r=0). Default is twice the mesh size of the foil
                    winding.
                mesh_size_outer : float
                    The mesh size towards the boundary of the domain. Default is 50 time the mesh size of the foil
                    winding.
                dist_min : float
                    The distance from the foil winding until which the mesh size of the foil winding is still used.
                dist_max : float
                    The distance from the foil winding from which the `mesh_size_outer` is used.

            Between `dist_min` and `dist_max` the mesh size changes linear.
        """
        self.foil_winding.add_to_gmsh()
        dist = min(self.foil_winding.width, self.foil_winding.height)

        mesh_size_inner = kwargs.pop('mesh_size_inner', 2 * self.foil_winding.mesh_size)
        mesh_size_outer = kwargs.pop('mesh_size_outer', 50 * self.foil_winding.mesh_size)
        dist_min = kwargs.pop('dist_min', dist)
        dist_max = kwargs.pop('dist_max', 10 * dist)

        if self.radius_factor is None:
            height = self.height
            width = self.width
            height_mid = self.foil_winding.center[1]
            h = min(self.foil_winding.height, height / 2)
            points = [Point(0, height_mid - 0.5 * height, 0),
                      Point(width, height_mid - 0.5 * height, 0),
                      Point(width, height_mid + 0.5 * height, 0),
                      Point(0, height_mid + 0.5 * height, 0),
                      Point(0, height_mid + h, 0, mesh_size_inner),
                      Point(0, height_mid - h, 0, mesh_size_inner)]
            lines = [Line(points[k], points[k + 1]) for k in range(-1, len(points) - 1)]
            for line in lines:
                line.add_to_gmsh()

            self._outer_curves = lines[1:4]
            self._curves = lines
        else:
            z_mid = self.foil_winding.west[1]  # z coordinate of mid point of circle
            radius = self.radius

            center_point = Point(0, z_mid, 0)
            bot_point = Point(0, z_mid - radius, 0)
            top_point = Point(0, z_mid + radius, 0)
            h = min(self.foil_winding.height, radius / 2)
            p1 = Point(0, z_mid + h, 0, mesh_size_inner)
            p2 = Point(0, z_mid - h, 0, mesh_size_inner)

            p1.add_to_gmsh()
            p2.add_to_gmsh()

            center_point.add_to_gmsh()
            bot_point.add_to_gmsh()
            top_point.add_to_gmsh()
            arc = CircleArc(top_point, bot_point, center_point)
            lines = [Line(top_point, p1), Line(p1, p2), Line(p2, bot_point)]

            arc.add_to_gmsh()
            for line in lines:
                line.add_to_gmsh()

            self._outer_curves = [arc, ]
            self._curves = [arc, ] + lines
            center_point.remove()

        surrounding = Surface(self.curves, self.foil_winding.outer_lines)
        surrounding.add_to_gmsh()

        tag_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(tag_distance, "PointsList", [self.foil_winding.points[0][0].tag,
                                                                      self.foil_winding.points[1][0].tag,
                                                                      self.foil_winding.points[0][-1].tag,
                                                                      self.foil_winding.points[1][-1].tag])
        gmsh.model.mesh.field.setNumbers(tag_distance, "CurvesList", [self.foil_winding.lines[2][0].tag,
                                                                      self.foil_winding.lines[2][-1].tag])

        tag_threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(tag_threshold, "InField", tag_distance)
        gmsh.model.mesh.field.setNumber(tag_threshold, "SizeMin", self.foil_winding.mesh_size)
        gmsh.model.mesh.field.setNumber(tag_threshold, "SizeMax", mesh_size_outer)
        gmsh.model.mesh.field.setNumber(tag_threshold, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(tag_threshold, "DistMax", dist_max)

        gmsh.model.mesh.field.setAsBackgroundMesh(tag_threshold)

        self.tag_distance = tag_distance
        self.tag_threshold = tag_threshold

        self._surrounding = surrounding

    @property
    def outer_curves(self) -> List[Curve]:
        """List with curves on the outer boundary, i.e. **not** on the z-axis."""
        return self._outer_curves

    @property
    def curves(self) -> List[Curve]:
        """List of all outer curves. Includes those on the z-axis."""
        return self._curves

    @property
    def outer_curve_bottom(self) -> Curve:
        """The outer curve on the bottom side."""
        return self._outer_curves[0]

    @property
    def outer_curve_right(self) -> Curve:
        """The outer curve on the right side."""
        if len(self._outer_curves) == 1:
            return self._outer_curves[0]
        return self._outer_curves[1]

    @property
    def outer_curve_top(self) -> Curve:
        """The outer curve on the top side."""
        if len(self._outer_curves) == 1:
            return self._outer_curves[0]
        return self._outer_curves[2]

    @property
    def surrounding(self) -> List[Curve]:
        """Surface of the surrounding, i.e. between foil windings and outer boundary."""
        return self._surrounding
