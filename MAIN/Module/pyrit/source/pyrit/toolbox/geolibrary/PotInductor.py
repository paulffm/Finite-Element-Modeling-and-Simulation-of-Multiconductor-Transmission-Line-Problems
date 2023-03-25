# coding=utf-8
"""Model of a pot inductor

.. sectionauthor:: Bundschuh
"""

from dataclasses import dataclass, field
from typing import Union, Tuple, NoReturn, List

import numpy as np

from pyrit import get_logger
from pyrit.geometry import Geo, Surface
from . import FoilWinding, Yoke

logger = get_logger(__name__)


@dataclass
class PotInductor(Geo):
    """Model representing a pot inductor.

    In a pot inductor, the windings are in the inside of a hollow cylinder. This model is axisymmetric. The dimensions
    of the model are mainly defines by the dimensions of the `windings` and the `yoke` that both have to be submitted to
    the constructor. The main setting of this class is the `placement` of the windings inside the winding window.

    Examples
    --------
    First, a foil winding and a yoke are created. These are then passed to the pot inductor.

    >>> foil_winding = FoilWinding(width=0.5, height=1)
    >>> yoke = Yoke(width=1, height=2)
    >>> pot_inductor = PotInductor(windings=foil_winding, yoke=yoke, placement=[0.5,0.2])

    This can then be used with the geometry module

    >>> from pyrit.geometry.Geometry import Geometry
    >>> geo = Geometry('pot inductor', show_gui=True)
    >>> with geo.geometry():
    >>>     pot_inductor.add_to_gmsh()
    """

    windings: FoilWinding  #: A foil winding object
    yoke: Yoke  #: A yoke object

    name: str = ''  #: Name of the foil winding
    origin_x: float = 0  #: x-coordinate of the origin
    origin_y: float = 0  #: y-coordinate of the origin
    origin_z: float = 0  #: z-coordinate of the origin
    anchor: Union[Geo.Anchor, Tuple[float]] = Geo.Anchor.WEST  #: Anchor of the foil winding.
    #: The placement of the windings in the window. [0,0] is the south west and [1,1] the north east corner.
    #: Corner radii are ignored.
    placement: List[float] = field(default_factory=lambda: [0.5, 0.5])

    filling_surface: Surface = field(default=None, init=False)  #: The surface between windings and yoke

    def __post_init__(self):
        super().__init__(self.name, self.origin_x, self.origin_y, self.origin_z, self.anchor)

        # Yoke and self need to have the same origin and anchor
        self.yoke.set_origin(self.origin_x, self.origin_y, self.origin_z)
        self.yoke.anchor = self.anchor

        # Check placement
        if not 0 < self.placement[0] < 1:
            logger.warning("First entry of placement out of bounds.")
        if not 0 < self.placement[1] < 1:
            logger.warning("Second entry of placement out of bounds.")

    def scale(self, scale: float) -> NoReturn:
        self.windings.scale(scale)
        self.yoke.scale(scale)

    def scale_x(self, scale: float):
        """Scale the model in x direction.

        Parameters
        ----------
        scale : float
            The factor to scale in x direction.
        """
        self.windings.scale_x(scale)
        self.yoke.scale_x(scale)

    def scale_y(self, scale: float):
        """Scale the model in y direction.

        Parameters
        ----------
        scale : float
            The factor to scale in y direction.
        """
        self.windings.scale_y(scale)
        self.yoke.scale_y(scale)

    def get_winding_coordinates(self) -> np.ndarray:
        """Get the coordinates of the center anchor of the winding.

        The placement of the windings in the winding window is determined by `placement`. This method computes and
        returns the coordinates of the center anchor of the windings that follows from the `placement`.

        Returns
        -------
        coordinates : np.ndarray
        """
        coordinates_zero_zero = self.yoke.south_west

        coordinates_zero_zero += np.array([self.yoke.get_thickness('thickness_inner'),
                                           self.yoke.get_thickness('thickness_bot'), 0])
        coordinates_zero_zero += np.array([0.5 * self.windings.width, 0.5 * self.windings.height, 0])

        coordinates_one_one = self.yoke.north_east
        coordinates_one_one -= np.array([self.yoke.get_thickness('thickness_outer'),
                                         self.yoke.get_thickness('thickness_top'), 0])
        coordinates_one_one -= np.array([0.5 * self.windings.width, 0.5 * self.windings.height, 0])

        difference = coordinates_one_one - coordinates_zero_zero

        result = coordinates_zero_zero
        result[:2] = result[:2] + np.multiply(self.placement, difference[:2])
        return result

    def add_to_gmsh(self, **kwargs) -> NoReturn:
        windings_kwargs = kwargs.pop('windings_kwargs', {})
        yoke_kwargs = kwargs.pop('yoke_kwargs', {})
        self.windings.set_origin(*self.get_winding_coordinates())
        self.windings.anchor = self.Anchor.CENTER

        self.windings.add_to_gmsh(**windings_kwargs)
        self.yoke.add_to_gmsh(**yoke_kwargs)

        self.filling_surface = Surface(self.yoke.inner_lines, self.windings.outer_lines)
        self.filling_surface.add_to_gmsh()

    def determine_coordinates(self, anchor: Union[Geo.Anchor, Tuple[float]]) -> np.ndarray:
        return self.yoke.determine_coordinates(anchor)
