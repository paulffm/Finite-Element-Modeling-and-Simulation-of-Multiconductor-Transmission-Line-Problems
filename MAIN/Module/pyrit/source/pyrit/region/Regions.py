# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 2021

Contains the class Regions

.. sectionauthor:: menzenbach
"""
from typing import List

from . import Regi


class Regions:
    """Class for Regions."""

    def __init__(self, *regis) -> None:
        """Constructor of Regions."""
        self.__regis = {}
        self.add_regi(*regis)

    def __iter__(self):
        return self.__regis.values().__iter__()

    def add_regi(self, *regis: Regi) -> None:
        """
        Add a Regi to the dictionary.

        Parameters
        ----------
        regis : Regi
            Regi instances which should be added to the dictionary.
        """
        for regi in regis:
            self.__regis[regi.ID] = regi

    def delete_regi(self, ID: int):
        """
        Delete a Regi from the dictionary.

        Parameters
        ----------
        ID : int
            The ID/key of the Regi which should be deleted from the dictionary.
        """
        # if ID does not exist in dict, None is returned
        self.__regis.pop(ID, None)

    def get_regi(self, ID: int) -> Regi:
        """
        Returns the Regi of with the given ID.

        Parameters
        ----------
        ID : int
           The ID/key of the Regi that should be returned.

        Returns
        -------
        Regi
            The Regi with the given ID.
        """
        return self.__regis.get(ID)

    def get_keys(self) -> List[int]:
        """
        Retuns all keys of the dictionary of all Regis.

        Returns
        -------
        List[int]
            The keys of the dictionary.
        """
        return list(self.__regis)

    def find_regions_of_boundary_condition(self, ID: int) -> List[int]:
        """Finds all regions that have the given boundary condition.

        Parameters
        ----------
        ID : int
            The ID of a boundary condition

        Returns
        -------
        indices_regions : List[int]
            A list of the IDs of the regions that have the boundary condition.
        """
        indices_regions = []
        for key in self.get_keys():
            regi = self.get_regi(key)
            if regi.bc == ID:
                indices_regions.append(regi.ID)
        return indices_regions

    def find_regions_of_excitations(self, ID: int) -> List[int]:
        """Finds all regions that have the given excitation.

        Parameters
        ----------
        ID : int
            The ID of an excitation

        Returns
        -------
        indices_regions : List[int]
            A list of the IDs of the regions that have the excitation.
        """
        indices_regions = []
        for key in self.get_keys():
            regi = self.get_regi(key)
            if regi.exci == ID:
                indices_regions.append(regi.ID)
        return indices_regions

    # def get_elem2mat(self, msh: Mesh, mats: Materials, mat_property: MatProperty):
    #     raise NotImplementedError
    # Check parameters
    # if not mats.get_ids():
    #     raise Exception(
    #         "Region does not contain any material information.")
    # if msh.node is None:
    #     raise Exception("Region does not contain any mesh information.")
    #
    # # Construct element-wise material property vector/matrix
    # num_elem = msh.elem2node.shape[0]   # number of elements
    #
    # # determine number of components of material property
    # # check if one material is isotrop
    # if any(mats.get_material(mat_id).get_property(mat_property).is_isotrop == True for mat_id in mats.get_ids()):
    #     mat_comp = 2
    # else:
    #     mat_comp = 1
    #
    # elem2mat = np.zeros((num_elem, mat_comp))
    #
    # for regi_id in self.get_keys():
    #     regi = self.get_regi(regi_id)
    #     mat = mats.get_material(regi.mat)
    #     mat_prop = mat.get_property(mat_property)
    #     # ?????????????????????????????? Niobe: find(prb.msh.elem2regi == r);
    #     idx_elem = 4
    #     if mat_prop.is_isotrop:
    #         mat_value = mat_prop.value
    #         if len(mat_value) == 1:
    #             mat_value *= np.ones((len(idx_elem), 1))
    #         elem2mat[idx_elem, :] = mat_value * np.ones((1, mat_comp))
    #     else:
    #         for cmp in range(mat_comp):
    #             # matValue = mat.mtcomp{cmp}(1).(property);
    #             pass
