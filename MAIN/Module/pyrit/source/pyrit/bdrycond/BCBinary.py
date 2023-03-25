# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:01:09 2021

Contains the abstract class BCBinary
.. sectionauthor:: Koray and menzenbach
"""

from numpy import ndarray

from . import BC


class BCBinary(BC):
    r"""Class of a binary boundary condition.

    Notes
    -----
    Binary boundary conditions mean :math:`\phi|_{\Gamma_{\mathrm{rep}}}=\gamma\phi|_{\Gamma_{\mathrm{prim}}}`, for a
    primary boundary :math:`\Gamma_{\mathrm{prim}}` and a replica boundary :math:`\Gamma_{\mathrm{rep}}`.
    """

    def __init__(self, value: float, name: str = ''):
        r"""
        Constructor of BCBinary.

        Parameters
        ----------
        value : float
            The value of $\gamma$
        name : str
            A name for the boundary condition
        """
        if isinstance(value, int):
            value = float(value)
        if not isinstance(value, float):
            raise ValueError(f"The type of value is {type(value)} but has to be float.")
        super().__init__(value, name)
        self._primary_nodes = None
        self._replica_nodes = None
        self._primary_edges = None
        self._replica_edges = None

    @property
    def primary_nodes(self) -> ndarray:
        """
        Getter property of primary_nodes

        Returns
        -------
        ndarray
            primary_nodes
        """
        return self._primary_nodes

    @primary_nodes.setter
    def primary_nodes(self, nodes: ndarray):
        """
        Setter property of primary_nodes

        Parameters
        ----------
        nodes : ndarray
            value the primary_nodes should have in an (N,)-ndarray

        Raises
        ------
        Exception
            if shape does not fit (N,)
        """
        if len(nodes.shape) != 1:
            raise Exception("parameter nodes should be of shape (N,)")
        self._primary_nodes = nodes

    @property
    def replica_nodes(self) -> ndarray:
        """
        Getter property of replica_nodes

        Returns
        -------
        ndarray
            replica_nodes
        """
        return self._replica_nodes

    @replica_nodes.setter
    def replica_nodes(self, nodes: ndarray):
        """
        Setter property of replica_nodes

        Parameters
        ----------
        nodes : ndarray
            value the replica_nodes should have in an (N,)-ndarray

        Raises
        ------
        Exception
            if shape does not fit (N,)
        """
        if len(nodes.shape) != 1:
            raise Exception("parameter nodes should be of shape (N,)")
        self._replica_nodes = nodes

    @property
    def primary_edges(self) -> ndarray:
        """
        Getter property of primary_edges

        Returns
        -------
        ndarray
            primary_edges
        """
        return self._primary_edges

    @primary_edges.setter
    def primary_edges(self, edges: ndarray):
        """
        Setter property of primary_edges

        Parameters
        ----------
        nodes : ndarray
            value the primary_edges should have in an (N,)-ndarray

        Raises
        ------
        Exception
            if shape does not fit (N,)
        """
        if len(edges.shape) != 1:
            raise Exception("parameter edges should be of shape (N,)")
        self._primary_edges = edges

    @property
    def replica_edges(self) -> ndarray:
        """
        Getter property of replica_edges

        Returns
        -------
        ndarray
            replica_edges
        """
        return self._replica_edges

    @replica_edges.setter
    def replica_edges(self, edges: ndarray):
        """
        Setter property of replica_edges

        Parameters
        ----------
        nodes : ndarray
            value the replica_edges should have in an (N,)-ndarray

        Raises
        ------
        Exception
            if shape does not fit (N,)
        """
        if len(edges.shape) != 1:
            raise Exception("parameter edges should be of shape (N,)")
        self._replica_edges = edges
