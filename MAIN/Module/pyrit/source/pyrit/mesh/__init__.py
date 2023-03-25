# coding=utf-8
"""
========================
Mesh (:mod:`pyrit.mesh`)
========================

.. currentmodule:: pyrit.mesh

In this package, classes for different kinds of meshes are defined, e.g. a triangular mesh for two-dimensional
domains or a tetrahedral mesh for three-dimensional domains. These classes serve mainly for bookkeeping. They store
the coordinates of all nodes and the definition of the higher dimensional entities (edges, triangles and tetrahedra).

Actual mesh classes
-------------------

.. autosummary::
    :toctree: autosummary/

    TriMesh
    AxiMesh
    TetMesh


Abstract class
--------------

.. autosummary::
    :toctree: autosummary/

    Mesh

"""

from .Mesh import Mesh
from .TriMesh import TriMesh
from .AxiMesh import AxiMesh
from .TetMesh import TetMesh

__all__ = ['Mesh', 'TriMesh', 'AxiMesh', 'TetMesh']
