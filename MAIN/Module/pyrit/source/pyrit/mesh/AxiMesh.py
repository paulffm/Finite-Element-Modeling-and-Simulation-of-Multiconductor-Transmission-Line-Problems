# coding=utf-8
"""
File containing the class AxiMesh.

.. sectionauthor:: bundschuh
"""

import warnings
from typing import NoReturn, Tuple
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt

from . import TriMesh


class AxiMesh(TriMesh):
    r"""Class representing a 2D-axisymmetric mesh in :math:`\rho z`-plane."""

    def __init__(self, node: np.ndarray, elem2node: np.ndarray, nodes_on_axis: np.ndarray = None):
        """
        Constructor of AxiMesh.

        Parameters
        ----------
        node :
        elem2node :
        nodes_on_axis : np.ndarray, optional
            Determines which nodes are on the z- axis. If not given, the indices are computed internally with the use of
            tolerance_on_axis (See method __calc_nodes_on_axis).
        """
        self.tolerance_on_axis: float = 1e-10  # Tolerance for when a node is considered as on the axis
        self.__check_axisymmetric(node)
        super().__init__(node, elem2node)
        if nodes_on_axis is None:
            self._nodes_on_axis = self.__calc_nodes_on_axis()
        else:
            self._nodes_on_axis = nodes_on_axis
        self._axi_element_area = None
        self._edges_on_axis = None
        self.__node_transformed = None
        self.__edge_length_transformed = None

    def __check_axisymmetric(self, node: np.ndarray) -> NoReturn:
        """
        Checks if the nodes all have a positive radius.

        Parameters
        ----------
        node : np.ndarray
            (N,2) array with radius and z-component in the first and second column, respectively. N nodes

        Raises
        ------
        ValueError
            When at least one node has a negative radius

        Returns
        -------
        None
        """
        if any(node[:, 0] < -1 * self.tolerance_on_axis):
            raise ValueError(
                "At least one node has a negative radius. The mesh is not a mesh for axisymmetric problems.")
        if any(node[:, 0] < 0):
            warnings.warn("At least one node has a radius that is slightly negative. This might be a rounding error.")

    def __calc_nodes_on_axis(self):
        """Calculates the indices of nodes that lie on the z-axis."""
        return np.nonzero(self.node[:, 0] < self.tolerance_on_axis)[0]

    def __calc_edges_on_axis(self):
        """Calculated the indices of the edges that lie on the z-axis."""
        elements_at_axis = set()
        for node_on_axis in self.nodes_on_axis:
            elements_at_axis.update(self.node2elem[node_on_axis])
        edges = set()
        for elem in elements_at_axis:
            edges.update(self.elem2edge[elem])

        edges_on_axis = []
        for edge in edges:
            if self.edge2node[edge][0] in self.nodes_on_axis and self.edge2node[edge][1] in self.nodes_on_axis:
                edges_on_axis.append(edge)
        return np.array(edges_on_axis)

    @property
    def node_transformed(self) -> np.ndarray:
        """The coordinates of the transformed nodes (s,z).

        Returns
        -------
        node_transformed : np.ndarray
            (T,2) array
        """
        if self.__node_transformed is None:
            self.__node_transformed = self.calc_node_transformed(self.node)
        return self.__node_transformed

    @property
    def edge_length_transformed(self) -> np.ndarray:
        """Length of all edges in the transformed coordinates (s,z).

        Returns
        -------
        np.ndarray
            (E,) array.
        """
        if self.__edge_length_transformed is None:
            self.__edge_length_transformed = self.calc_edge_length(self.edge2node, self.node_transformed)
        return self.__edge_length_transformed

    @property
    def nodes_on_axis(self) -> np.ndarray:
        """Array with indices of the nodes on the z-axis."""
        return self._nodes_on_axis

    @property
    def edges_on_axis(self):
        """Array with edge indices that are on the z-axis."""
        if self._edges_on_axis is None:
            self._edges_on_axis = self.__calc_edges_on_axis()
        return self._edges_on_axis

    @property
    def axi_element_area(self) -> np.ndarray:
        """Array with the area of each element in a mesh for edgefunctions."""
        if self._axi_element_area is None:
            self._axi_element_area = self.__calc_axi_element_area()
        return self._axi_element_area

    @staticmethod
    def generate_axi_mesh(msh: TriMesh, nodes_on_axis: np.ndarray = None):
        """
        Generates a object of AxiMesh from a TriMesh.

        Parameters
        ----------
        msh : TriMesh
        nodes_on_axis : np.ndarray, optional
            Determines which nodes are on the z- axis. If not given, the indices are computed internally with the use of
            tolerance_on_axis (See method __calc_nodes_on_axis).

        Returns
        -------
        AxiMesh
        """
        axi_mesh = AxiMesh(msh.node, msh.elem2node, nodes_on_axis)
        axi_mesh.node2regi = msh.node2regi
        axi_mesh.edge2regi = msh.edge2regi
        axi_mesh.elem2regi = msh.elem2regi

        return axi_mesh

    def plot_aximesh(self) -> Tuple[Figure, Axes]:
        """Plot the specific mesh for edgefunctions."""
        fig, ax = plt.subplots()

        for edge in range(self.num_edge):
            node1, node2 = self.edge2node[edge]
            if self.node[node1, 0] == self.node[node2, 0]:  # Nodes have same radius
                node_x = self.node[[node1, node2], 0]
                node_y = self.node[[node1, node2], 1]
                ax.plot(node_x, node_y, color='black')
            else:
                s1, z1 = self.node[node1, 0] ** 2, self.node[node1, 1]
                s2, z2 = self.node[node2, 0] ** 2, self.node[node2, 1]
                alpha = (z2 - z1) / (s2 - s1)
                beta = z1 - s1 * alpha

                rr = np.linspace(self.node[node1, 0], self.node[node2, 0], 50)
                zz = alpha * rr ** 2 + beta

                ax.plot(rr, zz, color='black')
        return fig, ax

    def __calc_axi_element_area(self):
        """Calculates the area of every element for the mesh needed for the mesh for Edgefunctions."""
        area = np.empty(self.num_elem)
        for k in range(self.num_elem):
            nodes = self.elem2node[k]
            r = self.node[nodes, 0]
            z = self.node[nodes, 1]

            matrix_b = np.array([[r[1] ** 2 - r[0] ** 2, r[2] ** 2 - r[0] ** 2], [z[1] - z[0], z[2] - z[0]]])
            abs_det_b = np.abs(np.linalg.det(matrix_b))

            alpha, beta, gamma = matrix_b[0, 0], matrix_b[0, 1], r[0] ** 2
            alpha_minus_beta = alpha - beta
            alpha_plus_gamma = alpha + gamma
            beta_plus_gamma = beta + gamma

            if alpha == 0:
                area[k] = abs_det_b * sum((-sqrt(gamma) / beta,
                                           - (-2 / 3 * sqrt(beta_plus_gamma) - 2 / 3 * gamma * sqrt(
                                               beta_plus_gamma) / beta) / beta,
                                           - 2 * gamma ** 1.5 / (3 * beta ** 2)))
            elif beta == 0:
                area[k] = abs_det_b * sum(((-2 * sqrt(gamma) + 4 * gamma ** 1.5 / (3 * alpha)) / (2 * alpha),
                                           sum((2 * sqrt(alpha_plus_gamma),
                                                2 * gamma / alpha * sum((-1 * gamma / sqrt(alpha_plus_gamma),
                                                                         - sqrt(alpha_plus_gamma))),
                                                2 / alpha * (gamma ** 2 / sqrt(alpha_plus_gamma) + 2 * gamma * sqrt(
                                                    alpha_plus_gamma) - 1 / 3 * alpha_plus_gamma ** 1.5))) / (
                                                   2 * alpha)))  # noqa: E126
            elif alpha != beta:
                area[k] = abs_det_b * sum((2 / 3 * (gamma ** 1.5 - alpha_plus_gamma ** 1.5) / (alpha * beta), 2 * (
                    sum((alpha_plus_gamma * sqrt(alpha_plus_gamma), - beta_plus_gamma * sqrt(beta_plus_gamma)))) / (
                                                   3 * beta * alpha_minus_beta)))  # noqa: E126
            else:
                area[k] = abs_det_b * sum(
                    (2 / 3 * (gamma ** 1.5 + alpha_plus_gamma ** 1.5) / (alpha * beta), sqrt(beta_plus_gamma) / beta))
        return area

    @staticmethod
    def calc_node_transformed(node: np.ndarray) -> np.ndarray:
        """Calculates the coordinates of the transformed nodes.

        Parameters
        ----------
        node : np.ndarray
            Coordinates of the nodes in (r,z) coordinate system

        Returns
        -------
        node_transformed : np.ndarray
            Coordinates of the transformed nodes in (s,z) coordinate system
        """
        return np.c_[node[:, 0] ** 2, node[:, 1]]
