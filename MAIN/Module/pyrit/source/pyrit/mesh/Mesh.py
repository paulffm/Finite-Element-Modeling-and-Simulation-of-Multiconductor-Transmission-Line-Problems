#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Contains the abstract class Mesh.

Created on Wed Mar 30 2021.

.. sectionauthor:: christ
@editor: bundschuh
"""

from abc import ABC, abstractmethod
import functools
from typing import List, TYPE_CHECKING, Type
import numpy as np

if TYPE_CHECKING:
    from pyrit.region import Regions
    from pyrit.material import Materials, MatProperty


class Mesh(ABC):
    """Abstract class for mesh object.

    Notes
    -----
    build a geometry and controls an external meshing algorithm (gmsh).
    """

    # __slots__ = ("_node", "_elem2node", "node2regi", "elem2regi")

    @abstractmethod
    def __init__(self, node: np.ndarray, elem2node: np.ndarray):
        """

        Constructor for the abstract class Mesh.

        Returns
        -------
        None.

        """
        self.node = node
        self.elem2node = elem2node
        self.node2regi = None
        self.elem2regi = None

    @property
    def node(self) -> np.ndarray:
        """
        Setter for node.

        Returns
        -------
        np.ndarray
            The nodes.

        """
        return self._node

    @node.setter
    def node(self, node):
        """
        Setter for node.

        Parameters
        ----------
        node : np.ndarray
            The nodes.

        Raises
        ------
        ValueError
            When the type of node is not np.ndarray.

        Returns
        -------
        None.

        """
        if not isinstance(node, np.ndarray):
            raise ValueError(
                f"Variable node is not a ndarray but a {str(type(node))}.")
        self._node = node

    @property
    def elem2node(self) -> np.ndarray:
        """
        Getter for elem2node.

        Returns
        -------
        np.ndarray
            elem2node.

        """
        return self._elem2node

    @elem2node.setter
    def elem2node(self, elem2node):
        """
        Setter for elem2node

        Parameters
        ----------
        elem2node : np.ndarray
            The relation bewteen elements and nodes.

        Raises
        ------
        ValueError
            When the type elem2node is not np.ndarray.

        Returns
        -------
        None.

        """
        if not isinstance(elem2node, np.ndarray):
            raise ValueError(
                f"Variable elem2node is not a ndarray but a {str(type(elem2node))}.")
        self._elem2node = elem2node

    @property
    @abstractmethod
    def bound2regi(self):
        """Getter for boundary entities to regions"""

    @bound2regi.setter
    @abstractmethod
    def bound2regi(self, bound2regi):
        """Setter for boundary entities to regions"""

    @property
    @abstractmethod
    def bound2node(self):
        """Getter for boundary entities to nodes"""

    @property
    def num_node(self) -> int:
        """
        Getter for num_node.

        Returns
        -------
        int
            The number of nodes in the mesh.

        """
        return self._node.shape[0]

    @property
    def num_elem(self) -> int:
        """
        Getter for num_elem

        Returns
        -------
        int
            The number of elements in the mesh.

        """
        return self._elem2node.shape[0]

    @abstractmethod
    def find_elemidx(self, coordinates: np.ndarray) -> np.ndarray:
        """Return the index of the element a given coordinate is located in.

        In the unlikely event that a coordinate is exactly on a node or an
        edge and thus, multiple elements have same shares in the coordinate,
        the element with the lowest index is returned.

        Parameters
        ----------
        coordinates : array_like
            Matrix of coordinates (or single coordinate) for which to find the
            element indices. Each row of 'cd' represents one coordinate.

        Returns
        -------
        np.ndarray
            array of element indices

        """

    @abstractmethod
    def find_nodeidx(self, coordinates: np.ndarray) -> np.ndarray:
        """Return the index of the node closest to a given coordinate.

        Choses the closest node based on Euclidian distance. In the unlikely
        event that a coordinate has numerically exact distances to nodes, the
        node with the lowest index is returned.

        Parameters
        ----------
        coordinates : np.ndarray
            Matrix of coordinates (or single coordinate) for which to find the
            closest node(s). Each row of 'cd' represents one coordinate.

        Returns
        -------
        np.ndarray
            single index or array of node indices

        """

    @staticmethod
    def calc_edge_length(edge2node: np.ndarray, node: np.ndarray) -> np.ndarray:
        """
        Computes the length of every edge.

        Parameters
        ----------
        edge2node : np.ndarray
            (E,2) array. E elements.
        node : np.ndarray
            (N,K) array. N nodes. K are the number of coordinated of one node and
            depends on the mesh.

        Returns
        -------
        np.ndarray
            (E,) array.

        """
        first_nodes = node[edge2node[:, 0], :]
        second_nodes = node[edge2node[:, 1], :]
        return np.sqrt(np.sum((first_nodes - second_nodes) ** 2, axis=1))

    @staticmethod
    def calc_face_area(face2node: np.ndarray, node: np.ndarray) -> np.ndarray:
        """
        Computes the length of every edge.

        Parameters
        ----------
        face2node : np.ndarray
            (F,3) array. F faces.
        node : np.ndarray
            (N,K) array. N nodes. K is the number of coordinates of one node and
            depends on the mesh.

        Returns
        -------
        np.ndarray
            (F,) array.

        """
        if node.shape[1] == 2:
            vector1 = node[face2node[:, 1], :] - node[face2node[:, 0], :]
            vector2 = node[face2node[:, 2], :] - node[face2node[:, 0], :]
            return 0.5 * np.abs(np.cross(vector1, vector2))
        if node.shape[1] == 3:
            vector1 = node[face2node[:, 1], :] - node[face2node[:, 0], :]
            vector2 = node[face2node[:, 2], :] - node[face2node[:, 0], :]
            return 0.5 * np.sqrt(np.sum(np.cross(vector1, vector2) ** 2, axis=1))
        raise ValueError(f"node has {node.shape[1]} columns but must have 2 or 3")

    @staticmethod
    def calc_tet_volume(tet2node: np.ndarray, node: np.ndarray) -> np.ndarray:
        """
        Computes the volume of every tetrahedron

        Parameters
        ----------
        tet2node : np.ndarray
            (T,4) array. T tetradedra
        node : np.ndarray
            (N,3) array. N nodes

        Returns
        -------
        np.ndarray
            (T,) array.

        """
        vector1 = node[tet2node[:, 1], :] - node[tet2node[:, 0], :]
        vector2 = node[tet2node[:, 2], :] - node[tet2node[:, 0], :]
        vector3 = node[tet2node[:, 3], :] - node[tet2node[:, 0], :]

        return 1 / 6 * np.abs(np.sum(vector1 * np.cross(vector2, vector3), axis=1))

    @staticmethod
    def calc_node2elem(elem2node: np.ndarray, num_node: int) -> List[List[int]]:
        """
        Calculates the node-to-element relation.

        Parameters
        ----------
        elem2node : np.ndarray
            (T,3) array with 3 nodes per Element
        num_node : int
            The number of nodes

        Returns
        -------
        List[List[int]]
        """
        nodes = [[] for _ in range(num_node)]
        for k in range(elem2node.shape[0]):
            for n in elem2node[k, :]:
                nodes[n].append(k)

        return nodes

    def value_per_element(self, regions: 'Regions', materials: 'Materials', prop_class: Type['MatProperty']):
        """Return a vector with one value per element.

        Evaluates the attribute `elem2regi` and returns a vector with the value of the given material property of the
        regions' material.

        Parameters
        ----------
        regions : Regions
            A regions object.
        materials : Materials
            A materials object.
        prop_class : Type[MatProperty]
            A class of MatProperty.

        Returns
        -------
        value_per_element : ndarray
            A vector with one value per element.
        """

        @functools.lru_cache
        def fun(regi_id, idx_elem):
            mat = materials.get_material(regions.get_regi(regi_id).mat)
            prop = mat.get_property(prop_class)
            value = mat.get_property(prop_class).value

            if isinstance(value, (float, int)):
                return value

            if not prop.is_homogeneous:
                raise NotImplementedError()
            if not prop.is_isotrop:
                raise NotImplementedError()
            if "element" in prop.keyword_args:
                return value(element=idx_elem)

            return value()

        iterable = (fun(self.elem2regi[elem_idx], elem_idx) for elem_idx in range(0, np.size(self.elem2regi)))
        return np.fromiter(iterable, float)

    @abstractmethod
    def boundary_elements(self) -> np.ndarray:
        """Returns an array with the indices of the entities at the boundary"""
