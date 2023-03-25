# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:44:03 2021

Contains the class Materials

.. sectionauthor:: bundschuh
"""
from warnings import warn
from typing import Type, Any, Iterable

from . import Mat, MatProperty


class Materials:
    """Class that manages materials of a problem"""

    __slots__ = ['__mats']

    def __init__(self, *mats: Mat):
        """
        Constructor of Material

        Parameters
        ----------
        *mats : Mat
            A number of Mats, possibly empty.

        Returns
        -------
        None.

        """
        self.__mats = {}
        if len(mats) > 0:
            self.add_material(*mats)

    def __iter__(self):
        return self.__mats.values().__iter__()

    def add_material(self, *mats: Mat):
        """
        Adds a number of Mats to Material.

        Parameters
        ----------
        *mats : Mat
            A number of Mats.

        Returns
        -------
        None.

        """
        for mat in mats:
            self.__mats[mat.ID] = mat

    def delete_material(self, *IDs: int):
        """
        Deletes a number of Mats specified by thier IDs from Material

        Parameters
        ----------
        *ids : int
            A number of IDs.

        Returns
        -------
        None.

        """
        for mat_ID in IDs:
            del self.__mats[mat_ID]

    def get_material(self, mat_ID: int):
        """
        Returns a number of Materials specified by their IDs

        Parameters
        ----------
        id : int
            An IDs.

        Raises
        ------
        KeyError
            When the ID of the Mat is unknown.

        Returns
        -------
        Mat
            A List of the Mats.

        """
        return self.__mats[mat_ID]

    def get_material_id_from_name(self, name: str) -> int:
        """
        Get the ID of a specific Mat from its name.

        The ID of the first Mat with the given name is returned. If multiple Mat's with the same name are present,
        a warning is given.

        Parameters
        ----------
        name : str
                The name of the requested Mat

        Returns
        -------
        id: int
            The id of the mat with the given name.
        """
        ids = self.get_ids()
        mat_names = [self.get_material(id).name for id in ids]
        if mat_names.count(name) > 1:
            msg = (f"The requested material named {name} is listed more than once. Returning the ID of the first Mat. "
                   f"Be aware of the global counter ID. Have you added the same material more than once or reused its "
                   f"name?")
            warn(msg)
        return ids[mat_names.index(name)]

    def get_material_from_name(self, name: str) -> Mat:
        """
        Get a specific Mat from its name. The first Mat with the given name is returned.

        Parameters
        ----------
        name : str
                The name of the requested Mat

        Returns
        -------
        mat: Mat
            The mat with the given name.
        """
        return self.get_material(self.get_material_id_from_name(name))

    def get_ids(self):
        """
        Returns the IDs of all Mats saved in Material

        Returns
        -------
        List[int]
            A list with the IDs of all Mats saved in Material.

        """
        return list(self.__mats.keys())

    def is_linear(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if the material is linear or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is linear. If not given, all properties are
            considered.

        Returns
        -------
        linear : bool
        """
        if prop_class is not None:
            return all((mat.is_linear(prop_class) for mat in self.__mats.values() if
                        mat.get_property(prop_class=prop_class) is not None))
        return all((mat.is_linear(prop_class) for mat in self.__mats.values()))

    def is_homogeneous(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if the material is homogeneous or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is homogeneous. If not given, all properties are
            considered.

        Returns
        -------
        homogeneous : bool
        """
        if prop_class is not None:
            return all((mat.is_homogeneous(prop_class) for mat in self.__mats.values() if
                        mat.get_property(prop_class=prop_class) is not None))
        return all((mat.is_homogeneous(prop_class) for mat in self.__mats.values()))

    def is_isotropic(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if the material is isotropic or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is isotropic. If not given, all properties are
            considered.

        Returns
        -------
        isotropic : bool
        """
        if prop_class is not None:
            return all((mat.is_isotropic(prop_class) for mat in self.__mats.values() if
                        mat.get_property(prop_class=prop_class) is not None))
        return all((mat.is_isotropic(prop_class) for mat in self.__mats.values()))

    def is_hysteretic(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if the material is hysteretic or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is hysteretic. If not given, all properties are
            considered.

        Returns
        -------
        hysteretic : bool
        """
        if prop_class is not None:
            return all((mat.is_hysteretic(prop_class) for mat in self.__mats.values() if
                        mat.get_property(prop_class=prop_class) is not None))
        return all((mat.is_hysteretic(prop_class) for mat in self.__mats.values()))

    def is_time_dependent(self, prop_class: Type[MatProperty] = None) -> bool:
        """Returns if the material is time dependent or not.

        Parameters
        ----------
        prop_class : Type[MatProperty], optional
            If given, it is only checked if the property of this type is time dependent. If not given, all properties
            are considered.

        Returns
        -------
        hysteretic : bool
        """
        if prop_class is not None:
            return any((mat.is_time_dependent(prop_class) for mat in self.__mats.values() if
                        mat.get_property(prop_class=prop_class) is not None))
        return any((mat.is_time_dependent(prop_class) for mat in self.__mats.values()))

    def update(self, name: str, val: Any, *, materials: Iterable = None, prop_class: Type[MatProperty] = None):
        """Update a value in the materials.

        Parameters
        ----------
        name : str
            Name of the argument in the MatProperties to update.
        val : Any
            The new value that is written on the MatProperties
        materials : Iterable, optional
            Iterable of materials to update values in. It can be integers (ID) or strings (material name)
        prop_class : Type[MatProperty], optional
            If given, set the value only in the MatProperty of given type. If not given, set the value for all.
            Default is None
        """
        if materials is None:
            mats = self.__mats.values()
        else:
            mats = (self.get_material(mat) if isinstance(mat, int) else self.get_material_from_name(mat) for mat in
                    materials)
        for mat in mats:
            if prop_class is None:
                mat.update(name, val, prop_class=prop_class)
            elif mat.get_property(prop_class) is not None:
                mat.update(name, val, prop_class=prop_class)
