# coding=utf-8
"""Class Surface. Describe a surface that can be made up in gmsh."""

from typing import List
from itertools import chain
import numpy as np
import gmsh


from . import PrimitiveGeo, Curve


class Surface(PrimitiveGeo):
    """A class representing a surface."""

    __slots__ = ('__ext_curves', '__int_curves')

    def __init__(self, ext_curves: List[Curve], *int_curves: List[Curve], physical_tag: int = None):
        """
        Construct an instance of Surface.

        Parameters
        ----------
        ext_curves: List[Curve]
            External curves.
        *int_curves: List[Curve]
            Internal curves. For each hole in the surface a list of curves that define the hole.
        physical_tag: int, optional
            The physical tag of the surface. The default is None.
        """
        super().__init__(2, physical_tag)

        self.__ext_curves = ext_curves

        if len(int_curves) > 0:
            if isinstance(int_curves[0], list):
                self.__int_curves = list(int_curves)
            else:
                self.__int_curves = list(list(int_curves))
        else:
            self.__int_curves = []

    @property
    def ext_curves(self) -> List[Curve]:
        """
        Get ext_curves.

        Returns
        -------
        List[Curve]
            External curves of the surface.
        """
        return self.__ext_curves

    @ext_curves.setter
    def ext_curves(self, ext_curves: List[Curve]):
        """
        Set ext_curves.

        Parameters
        ----------
        ext_curves: List[Curve]
            External curves.
        """
        self.__ext_curves = ext_curves

    @property
    def int_curves(self) -> List[Curve]:
        """
        Get int_curves.

        Returns
        -------
        List[Curve]
            Internal curves of the surface.
        """
        if len(self.__int_curves) == 0:
            return []
        return list(chain(*self.__int_curves))

    @int_curves.setter
    def int_curves(self, int_curves: List[Curve]):
        """
        Set int_curves.

        Parameters
        ----------
        int_curves : List[Curve]
            List of internal curves.
        """
        self.__int_curves = list(int_curves)

    @property
    def curves(self) -> List[Curve]:
        """
        Return a list of all curves.

        Returns
        -------
        List[Curve]
            All curves of the surface
        """
        return self.ext_curves + self.int_curves if self.has_holes else self.ext_curves

    @property
    def num_curves(self) -> int:
        """
        The number of curves.

        Returns
        -------
        int
            Number of curves.

        """
        return len(self.curves)

    @property
    def has_holes(self) -> bool:
        """
        Get has_holes.

        Returns
        -------
        bool
            True if surface has internal curves.
        """
        return bool(self.int_curves)

    # Superclass methods:
    @property
    def dim(self) -> int:
        """
        Get the dimension of the surface.

        Returns
        -------
        int
            Dimension of a surface.
        """
        return self._dim

    @property
    def added(self) -> bool:
        """
        Check if the surface is already added to gmsh.

        Returns
        -------
        bool
            True if surface has already been added to gmsh.
        """
        return bool(self.tag)

    def add_to_gmsh(self):
        """
        Add the surface to gmsh.

        Raises
        ------
        Warning
            When an object of class curve is part of the surface.

        Returns
        -------
        None.
        """
        if self.added:
            raise Warning("Surface is already added.")

        # check if curves are added in gmsh and add them:
        for cv in self.curves:
            if not cv.added:
                if type(cv) is Curve:  # pylint: disable=unidiomatic-typecheck
                    raise Warning("An object of class curve can not be added. "
                                  "It must either be an object of a subclass or it must be created by gmsh itself.")
                cv.add_to_gmsh()

        # add a closed wire for internal and external curves:
        # gmsh throws an exception if wire is not closed
        def tag(curve):  # return the tag of a curve
            return curve.tag

        try:
            tags = [c.tag for c in self.ext_curves]  # list(map(tag, self.ext_curves))
            wire_tags = [gmsh.model.occ.addWire(tags, checkClosed=True)]
            if self.has_holes:
                for cvs in self.__int_curves:
                    wire_tags.append(gmsh.model.occ.addWire(list(map(tag, cvs)), checkClosed=True))
        except Exception as ex:
            raise Warning("Curves must form a closed loop!") from ex

        # add surface as a plane surface in gmsh:
        self.tag = gmsh.model.occ.addPlaneSurface(wire_tags)

    def remove(self, recursive=True):  # pylint: disable=arguments-differ
        """
        Remove the Surface from GMSH.

        Returns
        -------
        None.

        """
        if self.added:
            gmsh.model.occ.remove([(self._dim, self.tag)], recursive=recursive)
            if recursive:
                self.reset_tag()
            else:
                self.tag = None
        else:
            raise Warning("The surface has not been added and thus can not be removed")

    def translate(self, dx: float, dy: float, dz: float):
        """
        Translate the Surface.

        Parameters
        ----------
        dx : float
            Translation in x-direction
        dy : float
            Translation in y-direction
        dz : float
            Translation in z-direction

        Returns
        -------
        None.
        """
        # translating the points of each curve (if curves are translated, points would be translated twice)
        points = []
        for cv in self.curves:
            for pt in cv.points:
                if pt not in points:
                    points.append(pt)

        for pt in points:
            pt.translate(dx, dy, dz)

        # translating the surface in gmsh
        if self.added:
            gmsh.model.occ.translate([(self.dim, self.tag)], dx, dy, dz)

    def rotate(self, angle_x: float, angle_y: float, angle_z: float, rotation_origin: np.ndarray = np.array([0, 0, 0])):
        """
        Rotate the Surface.

        Note the order of the rotation: first x-axis, second y-axis, third around z-axis.

        Parameters
        ----------
        angle_x : float
            Rotation angle around x-axis in radian
        angle_y : float
            Rotation angle around y-axis in radian
        angle_z : float
            Rotation angle around z-axis in radian
        rotation_origin : np.ndarray, optional
            The origin around which is rotated. By default (0,0,0) is used.

        Returns
        -------
        None.
        """
        # rotating points of surface: (points would be translated twice if curves are translated)
        points = []
        for cv in self.curves:
            for pt in cv.points:
                if pt not in points:
                    points.append(pt)

        for pt in points:
            pt.rotate(angle_x, angle_y, angle_z, rotation_origin)

        # rotating surface in gmsh:
        if self.added:
            x, y, z = rotation_origin
            gmsh.model.occ.rotate([(self.dim, self.tag)], x, y, z, 1, 0, 0, angle_x)
            gmsh.model.occ.rotate([(self.dim, self.tag)], x, y, z, 0, 1, 0, angle_y)
            gmsh.model.occ.rotate([(self.dim, self.tag)], x, y, z, 0, 0, 1, angle_z)

    def reset_tag(self) -> None:
        """
        Sets its own tag and the tags of all contained entities to None.

        Returns
        -------
        None
        """
        self.tag = None

        # reset tag of every curve part of surface:
        for cv in self.curves:
            cv.reset_tag()
