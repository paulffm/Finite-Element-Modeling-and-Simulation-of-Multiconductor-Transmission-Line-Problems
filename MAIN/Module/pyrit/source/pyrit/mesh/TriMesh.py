# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:27:37 2021

File containing the class TriMesh

.. sectionauthor:: bundschuh
"""

from typing import List, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.interpolate import griddata
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.artist import Artist
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from . import Mesh


class TriMesh(Mesh):
    """Class representing a triangular mesh in Cartesian coordinates."""

    # __slots__ = ('__elem2edge', '__edge2node', 'edge2regi', '__edge_length', '__elem_area')

    def __init__(self, node: np.ndarray, elem2node: np.ndarray):
        """
        Constructor of TriMesh

        Parameters
        ----------
        node : np.ndarray
            (N,2) array. Every line contains the x- and y-coordinate of a node. N nodes
        elem2node : np.ndarray
            (T,3) array. Every line contains the indices of the nodes that define an
            element (triangle). T elments

        Raises
        ------
        ValueError
            When node or elem2node is not a np.ndarray.
            When node does not have 2 columns

        Returns
        -------
        None.

        """
        super().__init__(node, np.sort(elem2node, axis=1))
        if node.shape[1] != 2:
            raise ValueError('The nodes do not have 2 coordinates')
        self.edge2regi = None

        # initializing private fields
        self.__elem2edge = None
        self.__edge2node = None
        self.__edge_length = None
        self.__elem_area = None
        self.__node2elem = None
        self.__edge2elem = None

    @property
    def elem2edge(self) -> np.ndarray:
        """
        element-to-edge relation.

        Returns
        -------
        np.ndarray
            (T,3) array. Every line contains the indices of the edges surrounding the element

        """
        if self.__elem2edge is None:
            self.__elem2edge = self.calc_elem2edge(self.elem2node, self.edge2node, self.num_node)
        return self.__elem2edge

    @property
    def edge2node(self) -> np.ndarray:
        """
        edge-to-node relation.

        Returns
        -------
        np.ndarray
            (E,2) array. Every line contains the indices of the nodes defining the edge. E edges

        """
        if self.__edge2node is None:
            self.__edge2node = self.calc_edge2node(self.elem2node)
        return self.__edge2node

    @edge2node.setter
    def edge2node(self, edge2node):
        """
        Setter for edge2node

        Parameters
        ----------
        edge2node : np.ndarray
            The relation bewteen edges and nodes.

        Raises
        ------
        ValueError
            When the type edge2node is not np.ndarray.

        Returns
        -------
        None.

        """
        if not isinstance(edge2node, np.ndarray):
            raise ValueError(
                f"Variable elem2node is not a ndarray but a {str(type(edge2node))}.")
        self.__edge2node = edge2node

    @property
    def edge_length(self) -> np.ndarray:
        """
        Length of all edges.

        Returns
        -------
        np.ndarray
            (E,) array.

        """
        if self.__edge_length is None:
            self.__edge_length = self.calc_edge_length(self.edge2node, self.node)
        return self.__edge_length

    @property
    def elem_area(self) -> np.ndarray:
        """
        Area of all elements

        Returns
        -------
        np.ndarray
            (E,) array.

        """
        if self.__elem_area is None:
            self.__elem_area = self.calc_face_area(self.elem2node, self.node)
        return self.__elem_area

    @property
    def num_edge(self) -> int:
        """
        Number of edges.

        Returns
        -------
        int
            Number of edges.

        """
        return self.edge2node.shape[0]

    @property
    def edge2elem(self) -> np.ndarray:
        """
        Edge-to-element relation

        Returns
        -------
        np.ndarray
            (E,2) array. Each line contains the indices of the (one or two) elements that contain
            the edge. For those edges that have only one elements, the second number is -1.

        """
        if self.__edge2elem is None:
            self.__edge2elem = self.calc_edge2elem(
                self.elem2edge, self.num_edge)
        return self.__edge2elem

    @property
    def node2elem(self):
        """Returns the node-to-element relation"""
        if self.__node2elem is None:
            self.__node2elem = self.calc_node2elem(self.elem2node, self.num_node)
        return self.__node2elem

    @property
    def bound2regi(self):
        return self.edge2regi

    @bound2regi.setter
    def bound2regi(self, bound2regi):
        self.edge2regi = bound2regi

    @property
    def bound2node(self):
        return self.edge2node

    def prj_el2nd(self, val2elem: np.ndarray) -> np.ndarray:
        """
        Projection of a element-based entity to a node-based entity

        Parameters
        ----------
        val2elem : np.ndarray
            (E,) array containing one value per element.

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
            (1 / 3 * np.ones(3 * self.num_elem),
             (np.repeat(np.arange(self.num_elem), 3), np.reshape(self.elem2node, 3 * self.num_elem))))
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
            (1 / 3 * np.ones(3 * self.num_elem),
             (np.repeat(np.arange(self.num_elem), 3), np.reshape(self.elem2node, 3 * self.num_elem))))
        return mat @ val2node

    def plot(self, entity_types: List[int] = (2,), entity_indices: List[np.ndarray] = None,
             label_nodes: bool = False, label_edges: bool = False,
             label_elements: bool = False, fig_info: Tuple[Figure, Axes] = None) -> Tuple[Figure, Axes]:
        """
        Plot the mesh.

        Parameters
        ----------
        entity_types : List[int], optional
            The dimension of the entities to be plotted. Possible dimensions are 0, 1 or 2.
            The default is [2].
        entity_indices : List[np.ndarray], optional
            The indices of the entities to be plotted. If given, this list must have the same
            length as `entity_types`. The default is None. Then, all entities are plotted.
        label_nodes : bool, optional
            If True, the index of every plotted node is displayed. The default is False.
        label_edges : bool, optional
            If True, the index of every plotted edge is displayed. The default is False.
        label_elements : bool, optional
            If True, the index of every plotted element is displayed. The default is False.

        Raises
        ------
        ValueError
            If the entity_type is not 0, 1 or 2.

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes

        """
        if fig_info is None:
            fig, ax = plt.subplots()
        else:
            fig = fig_info[0]
            ax = fig_info[1]

        for ent_type, k in zip(entity_types, range(len(entity_types))):
            if ent_type not in [0, 1, 2]:
                raise ValueError(
                    f"Entity type {str(ent_type)} is not allowed. Only 0,1,2 are allowed")
            if ent_type == 0:
                if entity_indices is None:
                    entity_indices = [np.arange(self.num_node)]
                ax.scatter(self.node[entity_indices[k], 0],
                           self.node[entity_indices[k], 1], color='black')
                if label_nodes:
                    for m in range(len(entity_indices[k])):
                        ax.annotate(str(entity_indices[k][m]), xy=(
                            self.node[entity_indices[k][m], 0], self.node[entity_indices[k][m], 1]), color='blue')
            if ent_type == 1:
                if entity_indices is None:
                    entity_indices = [np.arange(self.num_edge)]
                for m in entity_indices[k]:
                    node_x = self.node[self.edge2node[m, [0, 1]], 0]
                    node_y = self.node[self.edge2node[m, [0, 1]], 1]
                    # node_z = self.node[self.edge2node[m, [0, 1]], 2]
                    ax.plot(node_x, node_y, color='black')
                    if label_edges:
                        ax.annotate(str(entity_indices[k][m]), xy=(
                            np.mean(node_x), np.mean(node_y)), color='green')
                if label_nodes:
                    ind_points = np.unique(
                        self.edge2node[entity_indices[k], :])
                    for m in ind_points:
                        ax.annotate(str(m), xy=(
                            self.node[m, 0], self.node[m, 1]), color='blue')
            if ent_type == 2:
                if entity_indices is None:
                    for m in range(self.num_edge):
                        node_x = self.node[self.edge2node[m, [0, 1]], 0]
                        node_y = self.node[self.edge2node[m, [0, 1]], 1]
                        # node_z = self.node[self.edge2node[m, [0, 1]], 2]
                        ax.plot(node_x, node_y, color='black')
                    entity_indices = [np.arange(self.num_elem)]
                else:
                    for m in entity_indices[k]:
                        for n in range(3):
                            edge = self.elem2edge[m, n]
                            n1 = self.edge2node[np.abs(edge), 0]
                            n2 = self.edge2node[np.abs(edge), 1]

                            x = [self.node[n1, 0], self.node[n2, 0]]
                            y = [self.node[n1, 1], self.node[n2, 1]]
                            ax.plot(x, y, color='black')

                if label_elements:
                    for m in entity_indices[k]:
                        nodes = self.elem2node[m, :]
                        x = 1 / 3 * np.sum(self.node[nodes, 0])
                        y = 1 / 3 * np.sum(self.node[nodes, 1])
                        ax.annotate(str(m), xy=(x, y), color='red')
                if label_edges:
                    edges = np.unique(self.elem2edge[entity_indices[k], :])
                    for edge in edges:
                        node_x = self.node[self.edge2node[np.abs(edge), [
                            0, 1]], 0]
                        node_y = self.node[self.edge2node[np.abs(edge), [
                            0, 1]], 1]
                        ax.annotate(str(np.abs(edge)), xy=(
                            np.mean(node_x), np.mean(node_y)), color='green')
                if label_nodes:
                    nodes = np.unique(self.elem2node[entity_indices[k], :])
                    for node in nodes:
                        ax.annotate(str(node), xy=(
                            self.node[node, 0], self.node[node, 1]), color='blue')

        plt.title('The mesh')
        return fig, ax

    def plot_scalar_field(self, field: np.ndarray, plot_3d: bool = False, x_label: str = 'x axis',
                          y_label: str = 'y axis', z_label: str = 'z axis', title: str = '',
                          equal_axis_ratio: bool = True, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """
        Plot of a scalar field on the mesh. Requires one value per node.

        Parameters
        ----------
        field : np.ndarray
            (K,) array. A array with one value per node in the mesh. (K is the number of nodes or elements)
        plot_3d : bool, optional
            Indicates of the plot is 2- or 3-dimensional. Default is False
        x_label : str, optional
            Label of the x axis. Default is 'x axis'
        y_label : str, optional
            Label of the y axis. Default is 'y axis'
        z_label : str, optional
            Label of the z axis. Default is 'z axis'
        title : str, optional
            Title of the plot. Default is ''
        equal_axis_ratio: bool, optional
            If True, the axis ratios of x axis and y axis are equal

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        artist : Artist
            The artist
        """
        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        colorbar = kwargs.pop('colorbar', True)
        mesh = kwargs.pop('mesh', True)
        just_plot = kwargs.pop('just_plot', False)
        if just_plot:
            colorbar = False
            mesh = False

        if len(field.shape) == 2:
            if field.shape[1] == 1:
                field = field.flatten()
            else:
                raise ValueError(f"field has not the required format. It is {field.shape=} but should be (N,).")
        elif len(field.shape) > 2:
            raise ValueError(f"field has not the required format. It is {field.shape=} but should be (N,).")
        flag_nodal: bool
        if field.shape[0] == self.num_node:
            flag_nodal = True
        elif field.shape[0] == self.num_elem:
            flag_nodal = False
        else:
            raise ValueError(f"field has length {field.shape[0]} but should be {self.num_node}, the number of nodes.")
        if fig is None:
            fig = plt.figure()
        if ax is None:
            if flag_nodal:
                if plot_3d:
                    ax = fig.add_subplot(projection='3d')
                else:
                    ax = fig.add_subplot()
            else:
                ax = fig.add_subplot()
        # noinspection PyTypeChecker
        tri = Triangulation(self.node[:, 0], self.node[:, 1], self.elem2node)
        if flag_nodal:
            if plot_3d:
                art = ax.plot_trisurf(tri, field, cmap='viridis')  # cmap=plt.cm.CMRmap)
                ax.set_zlabel(z_label)
            else:
                art = ax.tripcolor(tri, field, shading='gouraud')
        else:
            art = ax.tripcolor(tri, facecolors=field)
        if colorbar:
            fig.colorbar(art, shrink=0.5, aspect=5, label=z_label)
        if mesh:
            ax.triplot(tri, color='black', lw=0.1)
        if equal_axis_ratio and not plot_3d:
            ax.set_aspect('equal', 'box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plt.title(title)
        return fig, ax, art

    def plot_vector_field(self, field: np.ndarray, norm_arrows: bool = False, color_arrows: bool = True,
                          x_label: str = 'x axis', y_label: str = 'y axis', title: str = '',
                          arrow_description: str = '', equal_axis_ratio: bool = True, **kwargs) -> \
            Tuple[Figure, Axes, Artist]:
        """
        Plot of a vector field on the mesh. Requires one vector per element.

        Parameters
        ----------
        field : np.ndarray
            (T,2) array. x- and y-component of a vector for every element. (T elements)
        norm_arrows : bool, optional
            If *True*, the length of all arrows is normalized. Default is *False*
        color_arrows : bool, optional
            If *True*, the arrows are colored according to their length. Default is *True*
        x_label : str, optional
            Label of the x axis. Default is 'x axis'
        y_label : str, optional
            Label of the y axis. Default is 'y axis'
        title : str, optional
            Title of the plot. Default is ''
        arrow_description : str, optional
            A legend for the arrows. Default is ''
        equal_axis_ratio: bool, optional
            If True, the axis ratios of x axis and y axis are equal

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        artist : Artist
            The artist

        See Also
        --------
        source.shapefunction.NodalShapeFunction.NodalShapeFunction.gradient: For generating the vector field from a
        scalar field
        """
        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        just_plot = kwargs.pop('just_plot', False)
        colorbar = kwargs.pop('colorbar', True)
        mesh = kwargs.pop('mesh', True)
        quiverkey = kwargs.pop('quiverkey', True)
        if just_plot:
            colorbar = False
            mesh = False
            quiverkey = False

        # Check format of field
        if len(field.shape) != 2:
            raise ValueError(f"field is of dimension {len(field.shape)} but should be of dimension 2")
        if field.shape[0] != self.num_elem or field.shape[1] != 2:
            raise ValueError(f"field has the shape {field.shape} but should be of shape {[self.num_elem, 2]}")

        if norm_arrows:
            length = np.hypot(field[:, 0], field[:, 1])
            arrows_x = field[:, 0] / length
            arrows_y = field[:, 1] / length
        else:
            arrows_x = field[:, 0]
            arrows_y = field[:, 1]
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        # noinspection PyTypeChecker
        if mesh:
            tri = Triangulation(self.node[:, 0], self.node[:, 1], self.elem2node)
            ax.triplot(tri, color='black', lw=0.5)
        midpoints = 1 / 3 * np.sum(self.node[self.elem2node], axis=1)
        if color_arrows:
            color_spec = np.hypot(field[:, 0], field[:, 1])
            quiv = ax.quiver(midpoints[:, 0], midpoints[:, 1], arrows_x, arrows_y, color_spec, units='xy', pivot='mid')
            if colorbar:
                fig.colorbar(quiv)
            if quiverkey:
                ax.quiverkey(quiv, 1.1, 1.05, 1, arrow_description, labelpos='E', coordinates='axes')
        else:
            quiv = ax.quiver(midpoints[:, 0], midpoints[:, 1], arrows_x, arrows_y, units='xy', pivot='mid')
            if quiverkey:
                ax.quiverkey(quiv, 0.9, 1.05, 1, arrow_description, labelpos='E', coordinates='axes')
        if equal_axis_ratio:
            ax.set_aspect('equal', 'box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plt.title(title)
        return fig, ax, quiv

    def plot_streamlines(self, field: np.ndarray, x_label: str = 'x axis', y_label: str = 'y axis', title: str = '',
                         equal_axis_ratio: bool = True, **kwargs) -> Tuple[Figure, Axes, Artist]:
        """
        Make a streamline plot of a given field.

        Parameters
        ----------
        field : np.ndarray
            (E,2) array. Two components of the field for every element E.
        x_label : str, optional
            Label of the x axis. Default is 'x axis'
        y_label : str, optional
            Label of the y axis. Default is 'y axis'
        title : str, optional
            Title of the plot. Default is ''
        kwargs : Any
            Keyword arguments passed to ax.streamplot
        equal_axis_ratio: bool, optional
            If True, the axis ratios of x axis and y axis are equal


        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        artist : Artist
            The artist

        Raises
        ------
        ValueError: When field has not the correct format
        """
        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        just_plot = kwargs.pop('just_plot', False)
        mesh = kwargs.pop('mesh', True)
        if just_plot:
            mesh = False

        if len(field.shape) != 2:
            raise ValueError(f"field is of dimension {len(field.shape)} but should be of dimension 2")
        if field.shape[0] != self.num_elem or field.shape[1] != 2:
            raise ValueError(f"field has the shape {field.shape} but should be of shape {[self.num_elem, 2]}")

        x_i = np.linspace(np.min(self.node[:, 0]), np.max(self.node[:, 0]), 200)
        y_i = np.linspace(np.min(self.node[:, 1]), np.max(self.node[:, 1]), 200)
        x_grid, y_grid = np.meshgrid(x_i, y_i)
        midpoints = 1 / 3 * np.sum(self.node[self.elem2node], axis=1)
        x_directions = griddata(midpoints, field[:, 0], (x_grid, y_grid), method='linear')
        y_directions = griddata(midpoints, field[:, 1], (x_grid, y_grid), method='linear')
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
        if mesh:
            tri = Triangulation(self.node[:, 0], self.node[:, 1], self.elem2node)
            ax.triplot(tri, color='gray', lw=0.5)
        stream = ax.streamplot(x_grid, y_grid, x_directions, y_directions, **kwargs)
        if equal_axis_ratio:
            ax.set_aspect('equal', 'box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plt.title(title)
        return fig, ax, stream

    def plot_equilines(self, field: np.ndarray, num_lines: int = 20, x_label: str = 'x axis',
                       y_label: str = 'y axis', title: str = '', equal_axis_ratio: bool = True, **kwargs) -> \
            Tuple[Figure, Axes, Artist]:
        """Plot equilines on the mesh.

        Parameters
        ----------
        field : ndarray
            The field to plot. Its size needs to be the number of nodes
        num_lines : int, optional
            The number of equilines. Default is 20
        x_label : str, optional
            Label of the x axis. Default is 'x axis'
        y_label : str, optional
            Label of the y axis. Default is 'y axis'
        title : str, optional
            Title of the plot. Default is ''
        equal_axis_ratio: bool, optional
            If True, the axis ratios of x axis and y axis are equal
        kwargs : Any
            Keyword arguments passed to ax.tricontour.

        Returns
        -------
        fig: Figure
            An object of matplotlib.figure.Figure.
        ax: Axes
            An object of matplotlib.axes.Axes
        artist : Artist
            The artist

        Raises
        ------
        ValueError: When field has not the correct format
        """
        fig = kwargs.pop('fig', None)
        ax = kwargs.pop('ax', None)
        just_plot = kwargs.pop('just_plot', False)
        mesh = kwargs.pop('mesh', True)
        if just_plot:
            mesh = False

        if field.size != self.num_node:
            raise ValueError(f"The field needs to be one value per node, but hat size {field.size}.")
        field = field.flatten()

        tri = Triangulation(self.node[:, 0], self.node[:, 1], self.elem2node)

        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)

        tpc = ax.tricontour(tri, field, num_lines, **kwargs)
        if mesh:
            ax.triplot(tri, color='black', lw=0.1)
        if equal_axis_ratio:
            ax.set_aspect('equal', 'box')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plt.title(title)

        return fig, ax, tpc

    @staticmethod
    def calc_elem2edge(elem2node: np.ndarray, edge2node: np.ndarray, num_node: int) -> np.ndarray:
        """
        Computes the element-to-edge realtion

        Parameters
        ----------
        elem2node : np.ndarray
            (T,3) array. T elements
        edge2node : np.ndarray
            (E,2) array. E edges
        num_node : int
            Number of nodes.

        Raises
        ------
        Warning
            If there is an error in the mesh.

        Returns
        -------
        np.ndarray
            (T,3) array.

        """
        # Get for every point all edges that use this point
        mat = sparse.coo_array((np.arange(edge2node.shape[0]), (edge2node[:, 0], edge2node[:, 1]))).tocsr()

        elem2edge = np.vstack([mat[elem2node[:, 0], elem2node[:, 1]], mat[elem2node[:, 1], elem2node[:, 2]],
                               mat[elem2node[:, 0], elem2node[:, 2]]]).transpose()

        return np.sort(elem2edge, axis=1)

    @staticmethod
    def calc_edge2node(elem2node: np.ndarray) -> np.ndarray:
        """
        Comutes the edge-to-node relation.

        Parameters
        ----------
        elem2node : np.ndarray
            (T,3) array. T elements

        Returns
        -------
        np.ndarray
            (E,2) array.

        """
        all_edges = np.r_[elem2node[:, [0, 1]],
                          elem2node[:, [1, 2]], elem2node[:, [0, 2]]]
        return np.unique(np.sort(all_edges), axis=0)

    @staticmethod
    def calc_edge2elem(elem2edge: np.ndarray, num_edges: int) -> np.ndarray:
        """
        Computes the edge-to-element relation

        Parameters
        ----------
        elem2edge : np.ndarray
            (T,3) array. T elements.
        num_edges : int
            Number of edges.

        Returns
        -------
        out : np.ndarray
            (E,2) array.

        """
        out = -1 * np.ones((num_edges, 2), dtype=int)
        lines = np.repeat(np.arange(elem2edge.shape[0]), 3)
        arr = np.abs(elem2edge).flatten()
        idx = np.argsort(arr)
        arr = arr[idx]
        lines = lines[idx]
        cols = np.abs(np.diff(arr) - 1)

        out[0, 0] = lines[0]
        out[arr[1:], cols] = lines[1:]

        return out

    def __point_in_element(self, elem_idx: int, p: np.ndarray) -> bool:
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
        bool, tuple
            (n,m): The coefficients to the vectors (b-a) and (c-a)
            1: out at edge between (b-a)
            2: out at edge between (c-a)
            3: out at edge between (c-b)

        Notes
        -----
        The three nodes defining the element are denoted as a, b and c
        """
        nodes = self.elem2node[elem_idx, :]
        a = self.node[nodes[0]]
        b = self.node[nodes[1]]
        c = self.node[nodes[2]]

        ap = p - a
        ab = b - a
        ac = c - a

        mat = np.c_[ab, ac]

        nm = np.linalg.solve(mat, ap)

        n = nm[0]
        m = nm[1]

        if n < 0:
            return 2
        if m < 0:
            return 1
        if n + m > 1:
            return 3
        return (n, m)

    # pylint: disable=arguments-differ
    def find_elemidx(self, coordinates: np.ndarray, return_entity: bool = False) -> np.ndarray:
        """
        For every point given in coordinates, this method finds the index of the element inside which the point lies.

        Parameters
        ----------
        coordinates : np.ndarray
            (N,2) array. Each line contains the x- and y-coordinates of a point.
        return_entity : bool, optional
            If True, additional information are returned: If a point lies on an
            edge or a node, the indices of these entities are returned. The default is False.

        Raises
        ------
        Exception
            When an internal error occures or when a requested point could not be found in the mesh.

        Returns
        -------
        elem_idx: np.ndarray
            (N,) array. Indices of the elements
        edge_idx: np.ndarray
            (N,) array. Indices of the edges a point may lie on.
            Contains -1 if the point does not lie on an edge. (optional)
        node_idx: np.ndarray
            (N,) array. Indices of the nodes a point may lie on.
            Contains -1 if the point does not lie on a node. (optional)

        Notes
        -----
        info/Algorithmen/find_elemidx

        """
        # adjust the shape of coordinates, if they do not match
        if len(coordinates.shape) == 1:
            coordinates = np.array(coordinates, ndmin=2)

        # Preparations
        local_elem2edge = np.abs(self.elem2edge)
        elem_idx = np.zeros(coordinates.shape[0], dtype=int)
        if return_entity:
            edge_idx = -1 * np.ones(coordinates.shape[0], dtype=int)
            node_idx = -1 * np.ones(coordinates.shape[0], dtype=int)

        for p in range(coordinates.shape[0]):  # For every point
            coordinate = coordinates[p, :]
            elems = np.arange(0, self.num_elem)
            elem = elems[0]

            while True:  # Iterate through the mesh
                # Check if the point lies in the current element
                out = self.__point_in_element(elem, coordinate)

                # Edges of the current element
                edges = local_elem2edge[elem, :]
                nodes = self.elem2node[elem]  # Nodes of the current element

                # Indices of first, second and third edge
                edge_n = [0, 0, 0]
                for k in range(3):
                    e2n = self.edge2node[edges[k]]
                    if nodes[0] in e2n and nodes[1] in e2n:
                        edge_n[0] = edges[k]
                    if nodes[0] in e2n and nodes[2] in e2n:
                        edge_n[1] = edges[k]
                    if nodes[1] in e2n and nodes[2] in e2n:
                        edge_n[2] = edges[k]

                # the point lies in the current element
                if isinstance(out, tuple):
                    if return_entity:  # Check if the point lies on an edge or a node
                        n = out[0]
                        m = out[1]
                        if np.isclose(n, 0):
                            if np.isclose(m, 0):
                                node_idx[p] = nodes[0]
                            elif np.isclose(m, 1):
                                node_idx[p] = nodes[2]
                            else:
                                edge_idx[p] = edge_n[1]
                        if np.isclose(m, 0):
                            if np.isclose(n, 1):
                                node_idx[p] = nodes[1]
                            elif not np.isclose(n, 0):
                                edge_idx[p] = edge_n[0]
                        if np.isclose(n + m, 1) and not np.isclose(n, 0) and not np.isclose(m, 0):
                            edge_idx[p] = edge_n[2]
                    break

                # Find the element for the next step
                for i in range(1, len(edge_n) + 1):
                    elems_to_edge = self.edge2elem[edge_n[out - i]]
                    if elems_to_edge[0] == elem:
                        new_elem = elems_to_edge[1]
                    else:
                        new_elem = elems_to_edge[0]

                    if new_elem == -1:  # The edge is at the boundary
                        continue

                    if new_elem not in elems:  # element has already been visited
                        continue
                    break
                if new_elem == -1:
                    raise Exception('An internal error occured')

                # Delete the current element from the list of not visited elements
                elems = elems[elems != elem]
                if elems.size == 0:  # No element was found that the point lies in
                    raise Exception(
                        f'It seems like the point at the coordinates x={coordinate[0]} and y={coordinate[1]} does not '
                        f'lie in the mesh.')
                if new_elem in elems:
                    elem = new_elem
                else:  # The element has already been visited. The method runs in circles.
                    raise Exception(
                        f'An internal error occured. The element with index {elem} has already been visited.')

            # Save the found element index
            elem_idx[p] = elem
        if return_entity:
            return (elem_idx, edge_idx, node_idx)
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

    def calc_normal_vector(self, edge: int):
        """Calculates the outer normal vector for the given edge.

        Parameters
        ----------
        edge : int
            The edge number.

        Returns
        -------
        n: ndarray
            The outer normal vector of the edge.
        """
        first_node, second_node = self.edge2node[edge]
        delta_r = self.node[second_node, 0] - self.node[first_node, 0]
        delta_z = self.node[second_node, 1] - self.node[first_node, 1]
        out = np.array([delta_z, -delta_r])
        nodes_element = self.elem2node[self.edge2elem[edge][0]]
        third_node = np.setdiff1d(nodes_element, self.edge2node[edge])
        sign = -1 * np.sign(np.dot(self.node[third_node] - self.node[first_node], out))

        return sign / np.linalg.norm(out) * out

    def boundary_elements(self) -> np.ndarray:
        """Gives the indices of all edges at the boundary."""
        return np.nonzero(self.edge2elem[:, 1] == -1)[0]
