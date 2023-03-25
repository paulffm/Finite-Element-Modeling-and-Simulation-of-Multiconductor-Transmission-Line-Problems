# -*- coding: utf-8 -*-
"""
Created on Fri May 14 12:47:57 2021

.. sectionauthor:: bundschuh
"""

from typing import List, Tuple, Type

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from . import Mesh


class TetMesh(Mesh):
    """Class representing a tetrahedral mesh."""

    def __init__(self, node: np.ndarray, elem2node: np.ndarray):
        super().__init__(node, elem2node)
        if node.shape[1] != 3:
            raise ValueError('The nodes do not have 3 coordinates')

        self.face2regi = None
        self.edge2regi = None

        # initializing private fields
        self.__elem2edge = None
        self.__elem2face = None
        self.__face2node = None
        self.__face2edge = None
        self.__edge2node = None
        self.__face2elem = None
        self.__edge_length = None
        self.__face_area = None
        self.__elem_vol = None

        self.__node2elem = None

    # %% Properties

    @property
    def elem2edge(self) -> np.ndarray:
        """
        Element-to-edge relation

        Returns
        -------
        elem2edge : np.ndarray
            (T,6) array. T elements and for every element the indices of the 6 surrounding edges
        """
        if self.__elem2edge is None:
            self.__elem2edge = self.calc_elem2edge(self.elem2node, self.edge2node, self.num_node)
        return self.__elem2edge

    @property
    def elem2face(self) -> np.ndarray:
        """
        Element-to-face relation

        Returns
        -------
        elem2face : np.ndarray
            (T,4) array. T elements and for every element the indices of the 4 surrounding faces
        """
        if self.__elem2face is None:
            self.__elem2face = self.calc_elem2face(self.elem2node, self.face2node, self.num_node)
        return self.__elem2face

    @property
    def face2node(self) -> np.ndarray:
        """
        Face-to-node relation

        Returns
        -------
        face2node : np.ndarray
            (F,3) array. F faces and for every face the indices of the defining nodes
        """
        if self.__face2node is None:
            self.__face2node = self.calc_face2node(self.elem2node)
        return self.__face2node

    @property
    def face2edge(self) -> np.ndarray:
        """
        Face-to-edge relation

        Returns
        -------
        face2edge : np.ndarray
            (F,3) array. F faces and for every face the indices of the surrounding edges
        """
        if self.__face2edge is None:
            self.__face2edge = self.calc_face2edge(self.face2node, self.edge2node, self.num_node)
        return self.__face2edge

    @property
    def edge2node(self) -> np.ndarray:
        """
        Edge-to-node relation

        Returns
        -------
        edge2node : np.ndarray
            (E,2) array. E edges and for every edge the defining nodes
        """
        if self.__edge2node is None:
            self.__edge2node = self.calc_edge2node(self.elem2node)
        return self.__edge2node

    @property
    def face2elem(self) -> np.ndarray:
        """
        Face-to-element relation

        Returns
        -------
        face2elem : np.ndarray
            (F,2) array. F faces and 2 (or 1) element per face. If only one element per faces exists
            the second entry in this line is -1.
        """
        if self.__face2elem is None:
            self.__face2elem = self.calc_face2elem(self.elem2face, self.num_face)
        return self.__face2elem

    @property
    def edge_length(self) -> np.ndarray:
        """
        The length of every edge

        Returns
        -------
        edge_length : np.ndarray
            (E,) array. E edges and the length of every edge
        """
        if self.__edge_length is None:
            self.__edge_length = self.calc_edge_length(self.edge2node, self.node)
        return self.__edge_length

    @property
    def face_area(self) -> np.ndarray:
        """
        The area of every face

        Returns
        -------
        face_area : np.ndarray
            (F,) array. F faces and the area of every face
        """
        if self.__face_area is None:
            self.__face_area = self.calc_face_area(self.face2node, self.node)
        return self.__face_area

    @property
    def elem_vol(self) -> np.ndarray:
        """
        The volume of every element

        Returns
        -------
        elem_vol : np.ndarray
            (T,) array. T elements and the volume of every element
        """
        if self.__elem_vol is None:
            self.__elem_vol = self.calc_tet_volume(self.elem2node, self.node)
        return self.__elem_vol

    @property
    def num_edge(self) -> int:
        """
        The number of edges in the mesh.

        Returns
        -------
        num_edge : int
            Number of edges
        """
        return self.edge2node.shape[0]

    @property
    def num_face(self) -> int:
        """
        The number of faces in the mesh.

        Returns
        -------
        num_face: int
            Number of faces
        """
        return self.face2node.shape[0]

    @property
    def node2elem(self):
        """Returns the node-to-element relation"""
        if self.__node2elem is None:
            self.__node2elem = self.calc_node2elem(self.elem2node, self.num_node)
        return self.__node2elem

    @property
    def bound2regi(self):
        return self.face2regi

    @bound2regi.setter
    def bound2regi(self, bound2regi):
        self.face2regi = bound2regi

    @property
    def bound2node(self):
        return self.face2node

    # %% Methods

    def prj_el2nd(self, val2elem: np.ndarray) -> np.ndarray:
        """
        Projection of a element-based entity to a node-based entity

        Parameters
        ----------
        val2elem : np.ndarray
            (E,) array containing one value per element

        Raises
        ------
        ValueError
            If the size of the input parameter does not match.

        Returns
        -------
        np.ndarray
            (N,) array containing one value per node.

        Notes
        -----
        The projection is implemented with a least-squares solution.

        """
        if not val2elem.shape[0] == self.num_elem:
            raise ValueError(
                f'Input parameter has the wrong size. It is {val2elem.shape[0]} but should be {self.num_elem}.')
        mat = sparse.coo_matrix(
            (1 / 4 * np.ones(4 * self.num_elem),
             (np.repeat(np.arange(self.num_elem), 4), np.reshape(self.elem2node, 4 * self.num_elem))))
        # noinspection PyTypeChecker
        val2node = lsqr(mat, val2elem, atol=1e-9)
        return val2node[0]

    def prj_nd2el(self, val2node: np.ndarray) -> np.ndarray:
        """
        Projection of a node-based entity to a element-based entity

        Parameters
        ----------
        val2node : np.ndarray
            (N,) array containing one value per node.

        Raises
        ------
        ValueError
            If the size of the input parameter does not match.

        Returns
        -------
        np.ndarray
            (E,) array containing one value per element.

        Notes
        -----
        The value per element is the mean of the values at its defining nodes.

        """
        if not val2node.shape[0] == self.num_node:
            raise ValueError(
                f'Input parameter has the wrong size. It is {val2node.shape[0]} but should be {self.num_node}.')

        mat = sparse.coo_matrix(
            (1 / 4 * np.ones(4 * self.num_elem),
             (np.repeat(np.arange(self.num_elem), 4), np.reshape(self.elem2node, 4 * self.num_elem))))
        return mat @ val2node

    def __point_in_element(self, elem_idx: int, p: np.ndarray) -> int:
        """
        Checks if the given point p is inside the element with index elem_idx.

        Parameters
        ----------
        elem_idx : int
            Index of the element.
        p : np.ndarray
            Coordinates of the point.

        Returns
        -------
        int
            0: point p is inside the element
            1: point p lies in direction of face 1
            2: point p lies in direction of face 2
            3: point p lies in direction of face 3
            4: point p lies in direction of face 4

        Notes
        -----
        The four nodes defining the element are denoted as a, b, c and d

        See Also
        --------
        info/Algorithmen/find_elemidx for nomenclature
        """
        nodes = self.elem2node[elem_idx, :]
        a = self.node[nodes[0]]
        b = self.node[nodes[1]]
        c = self.node[nodes[2]]
        d = self.node[nodes[3]]

        ap = p - a
        ab = b - a
        ac = c - a
        ad = d - a

        mat = np.c_[ab, ac, ad]
        nmo = np.linalg.solve(mat, ap)
        n = nmo[0]
        m = nmo[1]
        o = nmo[2]

        # Find the face in which direction p lies
        if n < 0 and m < 0 and o < 0:
            if o <= m and o <= n:
                return 4
            if n <= m and n <= o:
                return 2
            if m <= n and m <= o:
                return 3
        if n > 0 > m and o < 0:  # only n positive
            if m >= o:
                return 4
            return 3
        if n < 0 < m and o < 0:  # only m positive
            if n >= o:
                return 4
            return 2
        if n < 0 < o and m < 0:  # only o positive
            if m >= n:
                return 2
            return 3
        if n < 0 < m and o > 0:  # only n negative
            return 2
        if n > 0 > m and o > 0:  # only m negative
            return 3
        if n > 0 > o and m > 0:  # only o negative
            return 4
        if n + m + o > 1:
            return 1
        return 0

    def find_elemidx(self, coordinates: np.ndarray) -> np.ndarray:
        """
        For every point given in coordinates, this method finds the index of the element inside which the point lies.

        Parameters
        ----------
        coordinates : np.ndarray
            (N,3) array. Each line contains the x-, y- and z-coordinates of a point

        Raises
        ------
        Exception
            When an internal error occures or when a requested point could not be found in the mesh.

        Returns
        -------
        elem_idx: np.ndarray
            (N,) array. Indices of the elements

        Notes
        -----
        info/Algorithmen/find_elemidx
        """
        if len(coordinates.shape) == 1:
            coordinates = np.array(coordinates, ndmin=2)

        # Preparations
        local_elem2face = self.elem2face
        elem_idx = np.zeros(coordinates.shape[0], dtype=int)
        for p in range(coordinates.shape[0]):  # For every point
            coordinate = coordinates[p, :]
            elems = np.arange(0, self.num_elem)
            elem = elems[0]
            while True:  # Iterate through the mesh
                # Check if the point lies in the current element
                out = self.__point_in_element(elem, coordinate)
                if out == 0:
                    break

                faces = local_elem2face[elem, :]  # Edges of the current element
                nodes = self.elem2node[elem]  # Nodes of the current element

                # Indices of first, second, third and fourth face
                face_n = [0, 0, 0, 0]
                for k in range(4):
                    f2n = self.face2node[faces[k]]
                    for kk in range(4):
                        if nodes[kk] not in f2n:
                            face_n[kk] = faces[k]

                # Find the element for the next step
                for f in [face_n[out - 1], *face_n]:
                    elems_to_face = self.face2elem[f]
                    if elems_to_face[0] == elem:
                        new_elem = elems_to_face[1]
                    else:
                        new_elem = elems_to_face[0]

                    if new_elem == -1:  # The edge is at the boundary
                        continue
                    break
                if new_elem == -1:
                    raise Exception('An internal error occured')

                # Delete the current element from the list of not visited elements
                elems = elems[elems != elem]
                if elems.size == 0:  # No element was found that the point lies in
                    raise Exception(
                        f'It seems like the point at the coordinates x={coordinate[0]} and y={coordinate[1]} '
                        f'does not lie in the mesh.')
                if new_elem in elems:
                    elem = new_elem
                else:  # The element has already been visited. The method runs in circles.
                    raise Exception(
                        f'An internal error occured. The element with index {elem} has already been visited.')

            # Save the found element index
            elem_idx[p] = elem
        return elem_idx

    def find_nodeidx(self, coordinates: np.ndarray) -> np.ndarray:
        if len(coordinates.shape) == 1:
            coordinates = np.array(coordinates, ndmin=2)
        out = np.zeros(coordinates.shape[0])
        for k in range(coordinates.shape[0]):
            norms = np.linalg.norm(
                self.node - np.kron(np.ones((self.num_node, 1)), coordinates[k, :]), axis=1)
            out[k] = np.argmin(norms)
        return out

    def plot(self, entity_types=None, entity_indices: List[np.ndarray] = None,
             label_nodes: bool = False, label_edges: bool = False, label_faces: bool = False,
             label_elements: bool = False) -> Tuple[Figure, Type[Axes]]:
        """
        Plot the mesh.

        Parameters
        ----------
        entity_types : List[int], optional
            The dimension of the entities to be plotted. Possible dimensions are 0, 1, 2 or 3.
            The default is [3].
        entity_indices : List[np.ndarray], optional
            The indices of the entities to be plotted. If given, this list must have the same
            length as `entity_types`. The default is None. Then, all entities are plotted.
        label_nodes : bool, optional
            If True, the index of every plotted node is displayed. The default is False.
        label_edges : bool, optional
            If True, the index of every plotted edge is displayed. The default is False.
        label_faces : bool, optional
            If True, the index of every plotted face is displayed. The default is False.
        label_elements : bool, optional
            If True, the index of every plotted element is displayed. The default is False.

        Raises
        ------
        ValueError
            If the entity_type is not 0, 1, 2 or 3.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object.
        ax : matplotlib.axes.Axes
            The axes object.

        """
        if entity_types is None:
            entity_types = [3]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for ent_type, k in zip(entity_types, range(len(entity_types))):
            if ent_type not in [0, 1, 2, 3]:
                raise ValueError(
                    f"Entity type {str(ent_type)} is not allowed. Only 0,1,2 are allowed")

            if ent_type == 0:  # Plot nodes
                if entity_indices is None:
                    entity_indices = [np.arange(self.num_node)]
                ax.scatter(self.node[entity_indices[k], 0],
                           self.node[entity_indices[k], 1],
                           self.node[entity_indices[k], 2], color='black')
                if label_nodes:
                    self.__label_nodes(ax, entity_indices[k])

            if ent_type == 1:  # Plot edges
                if entity_indices is None:
                    entity_indices = [np.arange(self.num_edge)]
                for m in entity_indices[k]:
                    node_x = self.node[self.edge2node[m, [0, 1]], 0]
                    node_y = self.node[self.edge2node[m, [0, 1]], 1]
                    node_z = self.node[self.edge2node[m, [0, 1]], 2]
                    ax.plot(node_x, node_y, node_z, color='black')
                if label_edges:
                    self.__label_edges(ax, entity_indices[k])
                if label_nodes:
                    self.__label_nodes(ax, np.unique(self.edge2node[entity_indices[k], :]))

            if ent_type == 2:  # Plot faces
                if entity_indices is None:
                    for m in range(self.num_edge):
                        node_x = self.node[self.edge2node[m, [0, 1]], 0]
                        node_y = self.node[self.edge2node[m, [0, 1]], 1]
                        node_z = self.node[self.edge2node[m, [0, 1]], 2]
                        ax.plot(node_x, node_y, node_z)
                    entity_indices = [np.arange(self.num_face)]
                else:
                    for m in entity_indices[k]:
                        for n in range(3):
                            edge = self.face2edge[m, n]
                            n1 = self.edge2node[np.abs(edge), 0]
                            n2 = self.edge2node[np.abs(edge), 1]

                            x = [self.node[n1, 0], self.node[n2, 0]]
                            y = [self.node[n1, 1], self.node[n2, 1]]
                            z = [self.node[n1, 2], self.node[n2, 2]]
                            ax.plot(x, y, z, color='black')
                if label_faces:
                    self.__label_faces(ax, entity_indices[k])
                if label_edges:
                    self.__label_edges(ax, np.unique(self.face2edge[entity_indices[k], :]))
                if label_nodes:
                    self.__label_nodes(ax, np.unique(self.face2node[entity_indices[k], :]))

            if ent_type == 3:  # Plot elements
                if entity_indices is None:
                    for m in range(self.num_edge):
                        node_x = self.node[self.edge2node[m, [0, 1]], 0]
                        node_y = self.node[self.edge2node[m, [0, 1]], 1]
                        node_z = self.node[self.edge2node[m, [0, 1]], 2]
                        ax.plot(node_x, node_y, node_z)
                    entity_indices = [np.arange(self.num_elem)]
                else:
                    for m in entity_indices[k]:
                        for n in range(6):
                            edge = self.elem2edge[m, n]
                            n1 = self.edge2node[edge, 0]
                            n2 = self.edge2node[edge, 1]

                            x = [self.node[n1, 0], self.node[n2, 0]]
                            y = [self.node[n1, 1], self.node[n2, 1]]
                            z = [self.node[n1, 2], self.node[n2, 2]]
                            ax.plot(x, y, z, color='black')
                if label_elements:
                    self.__label_elements(ax, entity_indices[k])
                if label_faces:
                    self.__label_faces(ax, np.unique(self.elem2face[entity_indices[k], :]))
                if label_edges:
                    self.__label_edges(ax, np.unique(self.elem2edge[entity_indices[k], :]))
                if label_nodes:
                    self.__label_nodes(ax, np.unique(self.elem2node[entity_indices[k], :]))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        return fig, ax

    def __label_nodes(self, ax, nodes: np.ndarray):
        """
        Annotate the index of every given node in the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes of a plot.
        nodes : np.ndarray
            (K,) array with the indices of all nodes to be annotated
        """
        for node in nodes:
            ax.text(self.node[node, 0], self.node[node, 1], self.node[node, 2], str(node), color='blue')

    def __label_edges(self, ax, edges: np.ndarray):
        """
        Annotate the index of every given edge in the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes of a plot.
        edges : np.ndarray
            (K,) array with the indices of all edges to be annotated
        """
        for edge in edges:
            node_x = self.node[self.edge2node[np.abs(edge), [0, 1]], 0]
            node_y = self.node[self.edge2node[np.abs(edge), [0, 1]], 1]
            node_z = self.node[self.edge2node[np.abs(edge), [0, 1]], 2]
            ax.text(np.mean(node_x), np.mean(node_y), np.mean(node_z), str(edge), color='green')

    def __label_faces(self, ax, faces: np.ndarray):
        """
        Annotate the index of every given face in the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes of a plot.
        faces : np.ndarray
            (K,) array with the indices of all faces to be annotated
        """
        for face in faces:
            nodes = self.face2node[face, :]
            x = 1 / 3 * np.sum(self.node[nodes, 0])
            y = 1 / 3 * np.sum(self.node[nodes, 1])
            z = 1 / 3 * np.sum(self.node[nodes, 2])
            ax.text(x, y, z, str(face), color='red')

    def __label_elements(self, ax, elements: np.ndarray):
        """
        Annotate the index of every given elements in the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes of a plot.
        elements : np.ndarray
            (K,) array with the indices of all elements to be annotated
        """
        for element in elements:
            nodes = self.elem2node[element, :]
            x = 1 / 4 * np.sum(self.node[nodes, 0])
            y = 1 / 4 * np.sum(self.node[nodes, 1])
            z = 1 / 4 * np.sum(self.node[nodes, 2])
            ax.text(x, y, z, str(element), color='magenta')

    # %% Static methods

    @staticmethod
    def calc_elem2edge(elem2node: np.ndarray, edge2node: np.ndarray, num_node: int) -> np.ndarray:
        """
        Computes the element-to-edge relation.

        Parameters
        ----------
        elem2node : np.ndarray
            (T,4) array. T elements and 4 nodes per element
        edge2node : np.ndarray
            (E,2) array. E edges and 2 nodes per edge
        num_node : int
            Number of nodes in the mesh

        Returns
        -------
        elem2edge : np.ndarray
            (T,6) array. T elements and 6 edges per element

        """
        edges_of_nodes = [np.where(edge2node == k)[0] for k in range(num_node)]
        elem2edge = np.zeros((elem2node.shape[0], 6), dtype=np.int32)
        for m in range(elem2node.shape[0]):
            for n1, n2, k in zip([0, 0, 0, 1, 2, 3], [1, 2, 3, 2, 3, 1], [0, 1, 2, 3, 4, 5]):
                ind = np.intersect1d(edges_of_nodes[elem2node[m, n1]],
                                     edges_of_nodes[elem2node[m, n2]])
                elem2edge[m, k] = ind
        return np.sort(elem2edge, axis=1)

    @staticmethod
    def calc_elem2face(elem2node: np.ndarray, face2node: np.ndarray, num_node: int) -> np.ndarray:
        """
        Computes the element-to-face relation.

        Parameters
        ----------
        elem2node : np.ndarray
            (T,4) array. T elements and 4 nodes per element
        face2node : np.ndarray
            (F,3) array. F faces and 3 nodes per face
        num_node : int
            Number of nodes in the mesh.

        Returns
        -------
        elem2face : np.ndarray
            (T,4) array. T elements and 4 faces per element

        """
        # for every node the indices of all faces that use this node
        faces_of_nodes = [np.where(face2node == k)[0] for k in range(num_node)]
        elem2face = np.zeros((elem2node.shape[0], 4), dtype=np.int32)
        for m in range(elem2node.shape[0]):
            for n1, n2, n3 in zip([0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1]):
                ind = np.intersect1d(np.intersect1d(faces_of_nodes[elem2node[m, n1]],
                                                    faces_of_nodes[elem2node[m, n2]]), faces_of_nodes[elem2node[m, n3]])
                elem2face[m, n1] = ind
        return np.sort(elem2face, axis=1)

    @staticmethod
    def calc_face2node(elem2node: np.ndarray) -> np.ndarray:
        """
        Computes the face-to-node relation

        Parameters
        ----------
        elem2node : np.ndarray
            (T,4) array. T elements and 4 nodes per element

        Returns
        -------
        face2node : np.ndarray
            (F,3) array. F faces and 3 nodes per face

        """
        all_faces = np.r_[elem2node[:, [0, 1, 2]], elem2node[:, [0, 1, 3]],
                          elem2node[:, [0, 2, 3]], elem2node[:, [1, 2, 3]]]
        return np.unique(np.sort(all_faces), axis=0)

    @staticmethod
    def calc_face2edge(face2node: np.ndarray, edge2node: np.ndarray, num_node: int) -> np.ndarray:
        """
        Computes the face-to-edge relation.

        Parameters
        ----------
        face2node : np.ndarray
            (F,3) array. F faces and 3 nodes per face
        edge2node : np.ndarray
            (E,2) array. E edges and 2 nodes per edge
        num_node : int
            Number of nodes in the mesh

        Returns
        -------
        face2edge : np.ndarray
            (F,3) array. F faces and 3 edges per face

        """
        edges_of_nodes = [np.where(edge2node == k)[0] for k in range(num_node)]
        face2edge = np.zeros((face2node.shape[0], 3), dtype=np.int32)
        for m in range(face2node.shape[0]):
            for n1, n2 in zip([0, 1, 2], [1, 2, 0]):
                ind = np.intersect1d(edges_of_nodes[face2node[m, n1]],
                                     edges_of_nodes[face2node[m, n2]])
                face2edge[m, n1] = ind
        return np.sort(face2edge, axis=1)

    @staticmethod
    def calc_edge2node(elem2node: np.ndarray) -> np.ndarray:
        """
        Computes the edge-to-node relation

        Parameters
        ----------
        elem2node : np.ndarray
            (T,4) array. T elements and 4 nodes per element

        Returns
        -------
        edge2node : np.ndarray
            (E,2) array. E edges and 2 nodes per edge

        """
        all_edges = np.r_[elem2node[:, [0, 1]], elem2node[:, [0, 2]], elem2node[:, [0, 3]],
                          elem2node[:, [1, 2]], elem2node[:, [2, 3]], elem2node[:, [3, 1]]]
        return np.unique(np.sort(all_edges), axis=0)

    @staticmethod
    def calc_face2elem(elem2face: np.ndarray, num_face: int) -> np.ndarray:
        """
        Computes the face-to-elem relation.

        Parameters
        ----------
        elem2face : np.ndarray
            (T,4) array. T elements and 4 faces per element
        num_face : int
            Number of faces in the mesh

        Returns
        -------
        face2elem : np.ndarray
            (F,2) array. F faces and 2 (or 1) element per face. If only one element per faces exists
            the second entry in this line is -1.
        """
        out = -1 * np.ones((num_face, 2), dtype=int)
        for k in range(elem2face.shape[0]):
            for m in range(4):
                if out[elem2face[k, m], 0] == -1:
                    out[elem2face[k, m], 0] = k
                else:
                    out[elem2face[k, m], 1] = k
        return out

    def boundary_elements(self) -> np.ndarray:
        """Gives the indices of all faces at the boundary."""
        return np.nonzero(self.face2elem[:, 1] == -1)[0]
