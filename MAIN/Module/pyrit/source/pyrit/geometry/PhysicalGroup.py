# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

Contains the class Physical_Group

.. sectionauthor:: menzenbach
"""

import gmsh


class PhysicalGroup:
    """Class for physical group object that can be added directly to GMSH."""

    __slots__ = ('__tag', '__dim', 'name', 'entity_tags')

    def __init__(self, tag: int, dim: int, name: str = None):
        """
        Constructor of the class Physical_Group.

        Parameters
        ----------
        tag : int
            Unique tag of the physical group. If tag < 0 gmsh sets the tag, otherwise the tag should be > 0.
        dim : int
            Dimension of the physical group.
        name : str, optional
            Unique name of the physical group, by default None.

        Returns
        -------
        None.
        """
        self.tag = tag
        self.__dim = dim
        self.name = name
        self.entity_tags = set()

    @property
    def tag(self):
        """
        Get the unique tag of the physical group.

        Returns
        -------
        int
            Unique tag of the physical group.
        """
        return self.__tag

    @tag.setter
    def tag(self, tag: int):
        """
        Sets the tag with the given tag except the given tag is 0. Then an Exception is raised.

        Parameters
        ----------
        tag : int
            Unique tag for each physical group.

        Returns
        -------
        None.

        Raises
        ------
        ValueError
            ValueError if the given tag equals 0.
        """
        if tag == 0:
            raise ValueError(
                "The tag equals 0. Please change the tag to a number > 0 or -1 if gmsh should set the tag.")
        self.__tag = tag

    @property
    def dim(self) -> int:
        """
        Get the dimension of the Physical_Group object

        The following dimensions exist:
            0: Point
            1: Line
            2: Surface
            3: Volume

        Returns
        -------
        int
            The dimension of the physical group.
        """
        return self.__dim

    def add_entity(self, *entities):
        """
        Add a number of tags to entity_tags.

        Parameters
        ----------
        *entities : PrimitiveGeo
            a number of entities

        Returns
        -------
        None.

        Raises
        ------
        Exception
            Exception if the dimension of an entity does not match the dimension of the physical group.
        """
        for e in entities:
            if isinstance(e, int):
                self.entity_tags.add(e)
            elif e.dim != self.__dim:
                raise Exception(
                    f"Entity {str(e)} has dimension {str(e.dim)} but the physical group has dimension {self.dim}.")
            else:
                self.entity_tags.add(e.tag)

    def remove_entity(self, *entities):
        """
        Removes entities from the physical group.

        Parameters
        ----------
        entities : PrimitiveGeo
            A number of entities

        Returns
        -------
        None
        """
        for e in entities:
            if isinstance(e, int):
                self.entity_tags.discard(e)
            elif e.dim != self.__dim:
                raise Exception(
                    f"Entity {str(e)} has dimension {str(e.dim)} but the physical group has dimension {self.dim}.")
            else:
                self.entity_tags.discard(e.tag)

    def add_to_gmsh(self):
        """
        Add the Physical_Group object in GMSH.

        Returns
        -------
        None.

        Raises
        ------
        Exception
            Exception if there are no entities added to the physical group.
        """
        if self.entity_tags:
            self.__tag = gmsh.model.addPhysicalGroup(
                self.__dim, list(self.entity_tags), self.__tag)

            # set name
            if self.name is None:
                self.name = 'PhysicalGroup_' + str(self.__tag)
            gmsh.model.setPhysicalName(self.__dim, self.__tag, self.name)
        else:
            raise Exception(
                f"No entities in physical group {self.name}.")
