# coding=utf-8
"""Implementation of abtract class Geo.

.. sectionauthor:: Bundschuh
"""

from abc import ABC, abstractmethod
from typing import NoReturn, Union, Tuple
from enum import Enum
from numpy import ndarray
import numpy as np


class NotAddedToGmshError(Exception):
    """Custom Exception

    This Exception is raised if a method cannot be done because the method `add_to_gmsh()` has not been called yet.
    """


# pylint: disable=invalid-name
class BasicAnchor(Enum):
    """Basic anchors for Geo.

    They determine on which point the origin is set. Possible choices are `SOUTH`, `EAST`, `NORTH`, `WEST`,
    `SOUTH_WEST`, `SOUTH_EAST`, `NORTH_WEST`, `NORTH_EAST` and `CENTER`.
    """

    SOUTH = (0.5, 0)
    EAST = (1, 0.5)
    NORTH = (0.5, 1)
    WEST = (0, 0.5)
    SOUTH_WEST = (0, 0)
    SOUTH_EAST = (1, 0)
    NORTH_WEST = (0, 1)
    NORTH_EAST = (1, 1)
    CENTER = (0.5, 0.5)


class Geo(ABC):
    """Abstract class Geo.

    Collection of implementations for gmsh. Every `Geo` has a name, an own origin and a list of other `Geo` objects.
    """

    Anchor = BasicAnchor

    # @abstractmethod
    def __init__(self, name: str, origin_x: float = 0, origin_y: float = 0, origin_z: float = 0,
                 anchor: Union[Anchor, Tuple[float]] = Anchor.SOUTH_WEST):
        """Constructor

        Parameters
        ----------
        name : str
            Name of the Geo
        origin_x : float, optional
            The x-coordinate of the origin. Default 0
        origin_y : float, optional
            The y-coordinate of the origin. Default 0
        origin_z : float, optional
            The z-coordinate of the origin. Default 0
        """
        self.name = name
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.origin_z = origin_z
        self.anchor = anchor

        self.geo_list = []

    @abstractmethod
    def scale(self, scale: float) -> NoReturn:
        """Scale in all directions by a factor.

        Parameters
        ----------
        scale: float
            The factor to scale.
        """

    @abstractmethod
    def add_to_gmsh(self) -> NoReturn:
        """Adds the Entities of the `Geo` to gmsh.

        When this method is called, gmsh needs to be initialized. So at least

        .. code-block:: python

            gmsh.initialize(sys.argv)

        has to be called beforehand. However, it is recommended to use the context manager
        :py:class:`~pyrit.geometry.Geometry.Geometry` from the class `Geometry`.
        """

    @abstractmethod
    def determine_coordinates(self, anchor: Union[Anchor, Tuple[float]]) -> np.ndarray:
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

    def translate(self, dx: float, dy: float, dz: float) -> NoReturn:
        """Translate the `Geo`.

        Parameters
        ----------
        dx : float
            Translation in x direction.
        dy : float
            Translation in y direction.
        dz : float
            Translation in z direction.
        """
        old_origin = self.get_origin()
        new_origin = old_origin + np.array([dx, dy, dz])
        self.set_origin(*new_origin)
        for geo in self.geo_list:
            geo.translate(dx, dy, dz)

    def set_origin(self, origin_x: float, origin_y: float, origin_z: float) -> NoReturn:
        """Set all coordinates of the origin at one.

        Parameters
        ----------
        origin_x : float
            The x-coordinate of the origin.
        origin_y : float
            The y-coordinate of the origin.
        origin_z : float
            The z-coordinate of the origin.
        """
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.origin_z = origin_z

    def get_origin(self) -> ndarray:
        """Returns an array with the x, y and z coordinate of the origin.

        Returns
        -------
        origin : ndarray
            (3,) array with x, y and z coordinate of the origin.
        """
        return np.array([self.origin_x, self.origin_y, self.origin_z], dtype=float)

    @property
    def south_west(self) -> np.ndarray:
        """Coordinates of the south-west anchor"""
        return self.determine_coordinates(self.Anchor.SOUTH_WEST)

    @property
    def south(self) -> np.ndarray:
        """Coordinates of the south anchor"""
        return self.determine_coordinates(self.Anchor.SOUTH)

    @property
    def south_east(self) -> np.ndarray:
        """Coordinates of the south-east anchor"""
        return self.determine_coordinates(self.Anchor.SOUTH_EAST)

    @property
    def east(self) -> np.ndarray:
        """Coordinates of the east anchor"""
        return self.determine_coordinates(self.Anchor.EAST)

    @property
    def north_east(self) -> np.ndarray:
        """Coordinates of the north_east anchor"""
        return self.determine_coordinates(self.Anchor.NORTH_EAST)

    @property
    def north(self) -> np.ndarray:
        """Coordinates of the north anchor"""
        return self.determine_coordinates(self.Anchor.NORTH)

    @property
    def north_west(self) -> np.ndarray:
        """Coordinates of the north-west anchor"""
        return self.determine_coordinates(self.Anchor.NORTH_WEST)

    @property
    def west(self) -> np.ndarray:
        """Coordinates of the west anchor"""
        return self.determine_coordinates(self.Anchor.WEST)

    @property
    def center(self) -> np.ndarray:
        """Coordinates of the center anchor"""
        return self.determine_coordinates(self.Anchor.CENTER)
