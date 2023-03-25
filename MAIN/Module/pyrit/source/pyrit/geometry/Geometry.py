# coding=utf-8
"""
Created on Wed Mar 10 12:31:32 2021

.. sectionauthor:: bundschuh
"""
# pylint: disable-all
# flake8: noqa

import sys
from contextlib import contextmanager
import gmsh
import numpy as np
from scipy import sparse
from numba import jit, njit, prange, vectorize
from typing import Union, List, NoReturn
from inspect import signature
import warnings

from pyrit.material import Mat
from pyrit.bdrycond import BC
from pyrit.excitation import Exci
from pyrit.mesh import Mesh
from pyrit.mesh import TriMesh
from pyrit.mesh import TetMesh
from pyrit.region import Regi
from pyrit.region import Regions

# Specific Primitives
from . import Point, Line, Circle, CircleArc, Ellipse, PhysicalGroup

from pyrit import get_logger

logger = get_logger(__name__)


# region Decorators
def transform_pg(func):
    r"""
    Decorator that finds the parameter 'physical_group' in the \*args or \**kwargs. If the given value is an object of
    type PhysicalGroup it is replaced by its tag

    Parameters
    ----------
    func : function

    Returns
    -------
    wrapper: function
        The function with right physical group
    """

    def wrapper(*args, **kwargs):
        if 'physical_group' in kwargs:
            if isinstance(kwargs['physical_group'], PhysicalGroup):
                kwargs['physical_group'] = kwargs['physical_group'].tag
        else:
            params = signature(func).parameters
            for key, k in zip(params.keys(), range(len(args))):
                if key == 'physical_group':
                    if isinstance(args[k], PhysicalGroup):
                        args = list(args)
                        args[k] = args[k].tag

        val = func(*args, **kwargs)
        return val

    return wrapper


# endregion


class Geometry:
    mesh_versions = [1, 2.2]  # List of supported mesh versions

    # Dict of supported 2D meshing algorithms
    algorithms_2D = {'Automatic': 2, 'Delaunay': 5, 'Frontal-Delaunay': 6, 'BAMG': 7,
                     'Frontal-Delaunay for Quads': 8, 'Packing of Parallelograms': 9,
                     'Initial Mesh Only': 3, 'MeshAdapt': 1}
    # Dict of supported 3D meshing algorithms
    algorithms_3D = {'Delaunay': 1, 'Frontal': 4,
                     'HXT': 10, 'MMG3D': 7, 'Initial Mesh Only': 3}

    def __init__(self, model_name: str, mesh_version: float = 2.2, mesh_size_factor: float = 1,
                 min_element_size: float = 0, max_element_size: float = 1e22,
                 algorithm_2D: str = 'Frontal-Delaunay', algorithm_3D: str = 'Delaunay',
                 show_gui: bool = False, show_gmsh_info: bool = False):
        # Options for the model
        self.model_name = model_name
        self.mesh_version = mesh_version
        self.mesh_size_factor = mesh_size_factor
        self.min_element_size = min_element_size
        self.max_element_size = max_element_size
        self.algorithm_2D = algorithm_2D
        self.algorithm_3D = algorithm_3D
        self.show_gui = show_gui
        self.show_gmsh_info = show_gmsh_info

        # Lists of entities
        self.primitives = set()
        # self.points = []
        # self.curves = []
        # self.surfaces = []
        # self.volumes = []
        # self.geos = []

        # Dict of physical groups and relations between them and materials,
        # boundary conditions and excitations
        self.physical_groups = {}
        self.material_to_physical_group = {}
        self.boundary_condition_to_physical_group = {}
        self.excitation_to_physical_group = {}

        self.__writing_gmsh = False  # Will be True, while information can be passed to gmsh
        self.__physical_groups_added = False

    # %% Properties

    @property
    def mesh_version(self):
        return self.__mesh_version

    @mesh_version.setter
    def mesh_version(self, mesh_version: float):
        if mesh_version in Geometry.mesh_versions:
            self.__mesh_version = mesh_version
        else:
            raise ValueError('Mesh Version not supported.')

    @property
    def mesh_size_factor(self):
        return self.__mesh_size_factor

    @mesh_size_factor.setter
    def mesh_size_factor(self, mesh_size_factor):
        if mesh_size_factor > 0:
            self.__mesh_size_factor = mesh_size_factor
        else:
            raise ValueError(
                'The mesh size factor needs to be a positive number')

    @property
    def min_element_size(self):
        return self.__min_element_size

    @min_element_size.setter
    def min_element_size(self, min_element_size):
        if min_element_size >= 0:
            if '_Geometry__max_element_size' in self.__dict__:
                if min_element_size > self.__max_element_size:
                    raise ValueError(
                        'The minimum element size cannot be greater than the maximum element size')
            self.__min_element_size = min_element_size
        else:
            raise ValueError(
                'The minimum element size needs to be a positive number')

    @property
    def max_element_size(self):
        return self.__max_element_size

    @max_element_size.setter
    def max_element_size(self, max_element_size):
        if max_element_size > 0:
            if '_Geometry__min_element_size' in self.__dict__:
                if max_element_size < self.__min_element_size:
                    raise ValueError(
                        'The maximum element size cannot be smaller than the maximum element size')
            self.__max_element_size = max_element_size
        else:
            raise ValueError(
                'The minimum element size needs to be a positive number')

    @property
    def algorithm_2D(self):
        return self.__algorithm_2D

    @algorithm_2D.setter
    def algorithm_2D(self, algorithm_2D):
        if algorithm_2D in Geometry.algorithms_2D.keys():
            self.__algorithm_2D = algorithm_2D
        else:
            raise ValueError(
                f"The given algorithm '{algorithm_2D}' is not supported.")

    @property
    def algorithm_3D(self):
        return self.__algorithm_3D

    @algorithm_3D.setter
    def algorithm_3D(self, algorithm_3D):
        if algorithm_3D in Geometry.algorithms_3D.keys():
            self.__algorithm_3D = algorithm_3D
        else:
            raise ValueError(
                f"The given algorithm '{algorithm_3D}' is not supported.")

    # %% Methods

    def __create_physical_groups(self):
        """
        Creates the physical groups that are saved in self.physical_groups in gmsh.

        Returns
        -------
        None.

        """
        if not self.__physical_groups_added:
            gmsh.model.occ.synchronize()
            for group in self.physical_groups.values():  # Add physical groups
                group.add_to_gmsh()
            gmsh.model.occ.synchronize()
        self.__physical_groups_added = True

    def __set_options(self):
        if self.__writing_gmsh:
            gmsh.option.setNumber('Mesh.MshFileVersion', self.mesh_version)
            gmsh.option.setNumber('Mesh.MeshSizeFactor', self.mesh_size_factor)
            gmsh.option.setNumber('Mesh.MeshSizeMin', self.min_element_size)
            gmsh.option.setNumber('Mesh.MeshSizeMax', self.max_element_size)
            gmsh.option.setNumber(
                'Mesh.Algorithm', Geometry.algorithms_2D[self.algorithm_2D])
            gmsh.option.setNumber('Mesh.Algorithm3D',
                                  Geometry.algorithms_3D[self.algorithm_3D])
            if self.show_gmsh_info:
                gmsh.option.setNumber("General.Terminal", 1)
            else:
                gmsh.option.setNumber("General.Terminal", 0)

    def __initialize_gmsh(self):
        """
        Initializes gmsh, adds a model and sets the options

        Returns
        -------
        None.

        """
        # Initializing
        try:
            gmsh.initialize(sys.argv)
        except Exception as e:
            logger.info('Failed to initialize. I try again.')
            gmsh.initialize(sys.argv)

        self.__writing_gmsh = True
        gmsh.model.add(self.model_name)

        # Passing the options to gmsh
        self.__set_options()

        # Reset tags
        # for primitive in self.primitives:
        #     primitive.tag = None

        # Reset physical groups
        for group in self.physical_groups.values():
            group.entity_tags = set()

        self.__physical_groups_added = False

    def __finalize_gmsh(self):
        """
        Finalizes gmsh. If the option is set, the gui is shown

        Returns
        -------
        None.

        """
        if self.show_gui:
            self.open_gui()
        gmsh.finalize()

        self.__writing_gmsh = False

    def __mesh(self, dim: int):
        """
        Tells gmsh to mesh the geometry for given dimension

        Parameters
        ----------
        dim : int
            The dimension of the mesh. Must be in [0,1,2,3].

        Raises
        ------
        ValueError
            When the dimension is not in [0,1,2,3].

        Returns
        -------
        None.

        """
        if dim not in [0, 1, 2, 3]:
            raise ValueError(
                f"Dimension is {dim} but must be one of these values: {[0, 1, 2, 3]}")
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dim=dim)

    def __enter__(self):
        """Enter for context manager."""
        self.__initialize_gmsh()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Exit for context manager"""
        gmsh.model.occ.synchronize()
        self.__finalize_gmsh()

    @contextmanager
    def geometry(self):
        """
        Context manager for using gmsh. In the beginning it initialized gmsh and at the end it finalizes gmsh.


        Yields
        ------
        gmsh :
            The gmsh module

        """
        self.__initialize_gmsh()
        try:
            yield self
            gmsh.model.occ.synchronize()
        finally:
            self.__finalize_gmsh()

    def create_mesh(self, dim: int = 3):
        """
        Creates a mesh of dimensions dim (0, 1, 2 or 3) and according to the settings of Geometry
        When the method is called inside the context manager, the physical groups are added to gmsh
        and then the mesh is created (in gmsh). It can be changed afterwards, because gmsh does not
        get finalized.
        When the method is called outside the context manager, gmsh is initialized first, then the
        entities and physical groups are added, the mesh is created and last gmsh is finalized. The
        mesh cannot be changed afterwards.

        Parameters
        ----------
        dim : int, optional
            Dimension of the mesh. Has to be 0, 1, 2 or 3. The default is 3.

        Returns
        -------
        None.

        """
        self.__set_options()
        if self.__writing_gmsh:  # Inside the context manager geometry()
            self.__create_physical_groups()
            self.__mesh(dim)
        else:
            warnings.warn("You used create_mesh() outside of the context manager geometry(). Thus, it has no effect.")

    def refine_mesh(self):
        """
        Refines the mesh inside gmsh. This method only does something inside the
        context manager geometry().

        Returns
        -------
        None.

        """
        if self.__writing_gmsh:  # Inside the context manager geometry()
            gmsh.model.mesh.refine()
        else:
            warnings.warn("You used refine_mesh() outside of the context manager geometry(). Thus, it has no effect.")

    def get_mesh(self, dim: int = 3):
        """
        Return the mesh object
        """
        if self.__writing_gmsh:
            if dim == 3:
                return self.GmshMeshReader.read_tet_mesh()
            elif dim == 2:
                return self.GmshMeshReader.read_trimesh()
                # mesh = TriMesh(node, elem2node)
                # mesh.node2regi = self.GmshMeshReader.get_physical_groups_per_element(
                #     0, node, None)
                # mesh.edge2regi = self.GmshMeshReader.get_physical_groups_per_element(
                #     1, node, mesh.edge2node)
                # mesh.elem2regi = self.GmshMeshReader.get_physical_groups_per_element(
                #     2, node, elem2node)
                # return mesh
            elif dim == 1:
                return self.GmshMeshReader.read_linemesh()
            else:
                raise ValueError(f"Dimension is {dim} but must be one of these values: {[1, 2, 3]}")
        else:
            warnings.warn("You used get_mesh() outside of the context manager geometry(). Thus, it has no effect.")

    def get_regions(self):
        """
        Returns the regions.

        Returns
        -------
        None.

        """
        regis = dict()
        for key in self.physical_groups.keys():
            regis[key] = Regi(key, self.physical_groups[key].dim)

        for mat_id in self.material_to_physical_group.keys():
            for tag in self.material_to_physical_group[mat_id]:
                if regis[tag].mat is not None:
                    warnings.warn(f"Overriding the material of Regi with id={str(tag)}. " +
                                  f"Former material ID was {str(regis[tag].mat)}. New material ID is {str(mat_id)}.")
                regis[tag].mat = mat_id

        for bc_id in self.boundary_condition_to_physical_group.keys():
            for tag in self.boundary_condition_to_physical_group[bc_id]:
                if regis[tag].bc is not None:
                    warnings.warn(f"Overriding the bc of Regi with id={str(tag)}. " +
                                  f"Former bc ID was {str(regis[tag].bc)}. New bc ID is {str(bc_id)}.")
                regis[tag].bc = bc_id

        for exci_id in self.excitation_to_physical_group.keys():
            for tag in self.excitation_to_physical_group[exci_id]:
                if regis[tag].exci is not None:
                    warnings.warn(f"Overriding the exci of Regi with id={str(tag)}. " +
                                  f"Former exci ID was {str(regis[tag].exci)}. New exci ID is {str(exci_id)}.")
                regis[tag].exci = exci_id

        return Regions(*regis.values())

    def open_gui(self):
        """
        Opens the gmsh GUI to inspect the geometry (if possible, i.e. if inside the
        context manager geometry())

        Returns
        -------
        None.

        """
        try:
            if self.__writing_gmsh:  # Inside the context manager
                gmsh.model.occ.synchronize()
                gmsh.fltk.run()
            else:
                with self.geometry():
                    gmsh.fltk.run()
        except Exception as e:
            logger.error("Exception from gmsh: %s\nProceeding anyway", e.__str__())

    def save_geometry(self):
        pass

    def load_geometry(self):
        pass

    def add_physical_groups(self, *groups):
        for group in groups:
            if not isinstance(group, PhysicalGroup):
                raise TypeError(f"One given attribute has a not allowed type: {str(Geometry)}")
            self.physical_groups[group.tag] = group

    def create_physical_group(self, tag: int, dim: int, name: str = None):
        pg = PhysicalGroup(tag, dim, name)
        self.add_physical_groups(pg)
        return pg

    def add_material_to_physical_group(self, material: Union[int, Mat], *tags: Union[int, PhysicalGroup]) -> None:
        """
        Add a relation between a material and (possibly) multiple physical groups

        Parameters
        ----------
        material : Union[int, Mat]
            An instance of Mat or the ID of this Mat
        tags : Union[int, PhysicalGroup]
            A number of physical groups or tags of physical groups that have the specified material

        Raises
        ------
        TypeError
            If one input parameter has the wrong type

        Returns
        -------
        None
        """
        if isinstance(material, Mat):
            material = material.ID
        elif type(material) is not int:
            raise TypeError(f"Type of material is {str(type(material))} but only the types int and Mat are allowed.")
        if material not in self.material_to_physical_group.keys():
            self.material_to_physical_group[material] = list()

        for tag in tags:
            if isinstance(tag, PhysicalGroup):
                tag = tag.tag
            elif type(tag) is not int:
                raise TypeError(
                    f"Type of physical tag is {str(type(tag))} but only the types int and PhysicalGroup are allowed.")
            if tag not in self.physical_groups.keys():
                warnings.warn(f"There does not exist a physical group with tag = {str(tag)}.")
            self.material_to_physical_group[material].append(tag)

    def add_boundary_condition_to_physical_group(self, bc: Union[int, BC], *tags: Union[int, PhysicalGroup]) -> None:
        """
        Add a relation between boundary conditions and (possibly) multiple physical groups.

        Parameters
        ----------
        bc : Union[int, BC]
            An instance of BC or the ID of the instance
        tags : Union[int, PhysicalGroup]
            A number of physical groups or tags of physical groups that have the specified material

        Raises
        ------
        TypeError
            If one input parameter has the wrong type

        Returns
        -------
        None
        """
        if isinstance(bc, BC):
            bc = bc.ID
        elif type(bc) is not int:
            raise TypeError(f"Type of bc is {str(type(bc))} but only the types int and BC are allowed.")
        if bc not in self.boundary_condition_to_physical_group.keys():
            self.boundary_condition_to_physical_group[bc] = list()

        for tag in tags:
            if type(tag) is PhysicalGroup:
                tag = tag.tag
            elif type(tag) is not int:
                raise TypeError(
                    f"Type of physical tag is {str(type(tag))} but only the types int and PhysicalGroup are allowed.")
            if tag not in self.physical_groups.keys():
                warnings.warn(f"There does not exist a physical group with tag = {str(tag)}.")
            self.boundary_condition_to_physical_group[bc].append(tag)

    def add_excitation_to_physical_group(self, excitation: Union[int, Exci], *tags: Union[int, PhysicalGroup]) -> None:
        """
        Add a relation between excitations and (possibley) multiple physical groups.

        Parameters
        ----------
        excitation : Union[int Exci]
            An instance of Exci or the ID of the instance
        tags : Union[int, PhysicalGroup]
            A number of physical groups or tags of physical groups that have the specified material

        Raises
        ------
        TypeError
            If one input parameter has the wrong type

        Returns
        -------
        None
        """
        if isinstance(excitation, Exci):
            excitation = excitation.ID
        elif type(excitation) is not int:
            raise TypeError(
                f"Type of excitation is {str(type(excitation))} but only the types int and Exci are allowed.")
        if excitation not in self.excitation_to_physical_group.keys():
            self.excitation_to_physical_group[excitation] = list()
        for tag in tags:
            if type(tag) is PhysicalGroup:
                tag = tag.tag
            elif type(tag) is not int:
                raise TypeError(
                    f"Type of physical tag is {str(type(tag))} but only the types int and PhysicalGroup are allowed.")
            if tag not in self.physical_groups.keys():
                warnings.warn(f"There does not exist a physical group with tag = {str(tag)}.")
            self.excitation_to_physical_group[excitation].append(tag)

    def reset_tags(self):
        for primitive in self.primitives:
            primitive.reset_tag()

    def add_primitives(self, *primitives):
        r"""
        Adds instances of PrimitiveGeo to the Geometry.
        If called inside the context manager geometry() the instances are instantly added to gmsh.

        Parameters
        ----------
        primitives : PrimitiveGeo
            A number of instances of PrimitiveGeo.

        Raises
        ------
        Warning
            When an object in \*primitives is not an instance of PrimitiveGeo.

        Returns
        -------
        None.
        """
        if self.__writing_gmsh:  # Inside the context manager
            for primitive in primitives:
                if not primitive.added:
                    primitive.add_to_gmsh()
                if primitive.physical_tag is not None:
                    if primitive.physical_tag in self.physical_groups.keys():
                        self.physical_groups[primitive.physical_tag].add_entity(primitive)
                    else:
                        warnings.warn(
                            f"There does not exist a physical group with tag {primitive.physical_tag}." +
                            "Thus, the entity cannot be added to the physical group")
                self.primitives.add(primitive)
            gmsh.model.occ.synchronize()
        else:
            warnings.warn(
                "You used add_primitives() outside of the context manager geometry(). Thus, it has no effect.")

    def remove_primitives(self, *primitives):
        r"""
        Removes instances of PrimitiveGeo from Geometry.

        Parameters
        ----------
        primitives : PrimitiveGeo
            A number of instances of PrimitiveGeo.

        Raises
        ------
        Warning
            When an object in \*primitives is not an instance of PrimitiveGeo.

        Returns
        -------
        None.

        """
        for primitive in primitives:
            if self.__writing_gmsh:  # Inside the context manager
                if primitive.added:
                    primitive.remove()
            if primitive.physical_tag is not None:
                if primitive.physical_tag in self.physical_groups.keys():
                    try:
                        self.physical_groups[primitive.physical_tag].remove_entity(primitive)
                    except Exception as ex:
                        warnings.warn(f"PhysicalGroup raised an Exception:\n\t{str(ex)}\n" +
                                      "The Exception is not fatal. So it is proceeded")
            self.primitives.discard(primitive)
            primitive.reset_tag()

    def add_geos(self, *geos):
        pass

    def remove_geos(self, *geos):
        pass

    # %% Methods for directly adding primitiveGeo

    @transform_pg
    def add_point(self, x: float, y: float, z: float, mesh_size: float = None,
                  physical_group: Union[int, PhysicalGroup] = None):
        point = Point(x, y, z, mesh_size=mesh_size, physical_tag=physical_group)
        self.add_primitives(point)
        return point

    @transform_pg
    def add_line(self, point1: Point, point2: Point, physical_group: Union[int, PhysicalGroup] = None):
        line = Line(point1, point2, physical_group)
        self.add_primitives(line)
        return line

    class GmshMeshReader:

        @staticmethod
        def read_trimesh_2():
            node_tag, node_coords, _ = gmsh.model.mesh.getNodes()
            element_types, element_tags, node_tags_elements = gmsh.model.mesh.getElements()
            node_tag -= 1
            srt = np.argsort(node_tag)
            if not np.all(node_tag[srt] == np.arange(len(node_tag))):
                raise Exception(f"The node tags from gmsh are not as expected.")

            ind_points = np.where(element_types == 15)[0][0]  # Indices of points. 15 is the element type of points
            ind_edges = np.where(element_types == 1)[0][0]  # Indices of edges. 1 is the element type of edge
            ind_elements = np.where(element_types == 2)[0][0]  # Indices of triangles. 2 is the element type of triangle

            mesh_elem2node = node_tags_elements[ind_elements].reshape(-1, 3) - 1
            node = node_coords.reshape(-1, 3)
            mesh_node = node[node_tag[np.unique(mesh_elem2node)], 0:2]

            mesh = TriMesh(mesh_node, mesh_elem2node)
            return mesh

        @staticmethod
        def read_trimesh():
            # Read from gmsh
            node_tag, node, _ = gmsh.model.mesh.getNodes()
            element_types, element_tags, node_tags_elements = gmsh.model.mesh.getElements()
            if len(element_types) == 0:
                raise Exception("Cannot read mesh because no elements exist. Have you created the mesh?")

            # Getting nodes
            num_node = int(len(node) / 3)  # Number of nodes
            node = np.reshape(node, (num_node, 3))
            # Coordinates of nodes. x-coordinate in first column and y-coordinate in second column
            node = node[:, 0:2]
            node_tag = node_tag - np.ones(len(node_tag))
            node_tag = node_tag.astype('int')
            # node[node_tag, :] = node
            np.put_along_axis(node, np.c_[node_tag, node_tag], node, axis=0)

            # Getting tags
            ind_points = np.where(element_types == 15)[0][0]  # Indices of points. 15 is the element type of points
            ind_edges = np.where(element_types == 1)[0][0]  # Indices of edges. 1 is the element type of an edge
            ind_elements = np.where(element_types == 2)[0][0]  # Indices of triangles. 2 is the
            # element type of an triangle
            point_tags = element_tags[ind_points]
            edge_tags = element_tags[ind_edges]
            triangle_tags = element_tags[ind_elements]

            # Getting edges

            edges = np.array(node_tags_elements[ind_edges])
            num_edges = int(len(edges) / 2)  # Number of edges
            edge2node = np.reshape(edges, (num_edges, 2)) - 1
            # Connection between elements and nodes. Each line contains the indices of the contained nodes
            edge2node = edge2node.astype('int')
            edge2node = np.sort(edge2node, axis=1)

            # Getting elements (triangles)

            elements = np.array(node_tags_elements[ind_elements])
            num_elem = int(len(elements) / 3)  # Number of elements
            elem2node = np.reshape(elements, (num_elem, 3)) - np.ones((num_elem, 1))
            # Connection between elements and nodes. Each line contains the indices of the contained nodes
            elem2node = elem2node.astype('int')

            msh = TriMesh(node, elem2node)
            Geometry.GmshMeshReader.remove_unused(msh, TriMesh.calc_node2elem(elem2node, num_node))

            # Physical Groups of nodes
            pg = gmsh.model.getPhysicalGroups(dim=0)
            msh.node2regi = -1 * np.ones(num_node, dtype=int)
            for k in range(len(pg)):
                node_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind = node_tags_elements[ind_points][
                          np.array([np.where(point_tags == m)[0] for m in node_tags_pg]).flatten()] - 1
                msh.node2regi[ind] = pg[k][1]

            # Physical groups of edges
            pg = gmsh.model.getPhysicalGroups(dim=1)
            msh.edge2regi = -1 * np.ones(msh.num_edge, dtype=int)
            for k in range(len(pg)):
                edge_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind_tmp = np.array([np.where(edge_tags == m)[0] for m in edge_tags_pg]).flatten()
                tmp = edge2node[ind_tmp, :]
                ind = np.array([np.intersect1d(np.where(msh.edge2node[:, 0] == tmp[m, 0]),
                                               np.where(msh.edge2node[:, 1] == tmp[m, 1]))[0] for m in range(len(tmp))])
                msh.edge2regi[ind] = pg[k][1]

            # Physical groups of elements
            pg = gmsh.model.getPhysicalGroups(dim=2)
            msh.elem2regi = -1 * np.ones(num_elem, dtype=int)
            mat = sparse.coo_array(
                (np.arange(triangle_tags.size), (triangle_tags, np.zeros_like(triangle_tags)))).tocsr()
            for k in range(len(pg)):
                elem_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind = mat[elem_tags_pg, [0]]
                msh.elem2regi[ind] = pg[k][1]

            for node2elem in msh.node2elem:
                if not len(node2elem) > 0:
                    raise Exception("There are nodes in the mesh that don't belong to an element.")
            return msh

        @staticmethod
        def read_linemesh():
            # Read from gmsh
            node_tag, node, _ = gmsh.model.mesh.getNodes()
            element_types, element_tags, node_tags_elements = gmsh.model.mesh.getElements(
                dim=1)

            # Getting nodes
            num_node = int(len(node) / 3)  # Number of nodes
            node = np.reshape(node, (num_node, 3))
            # Coordinates of nodes. x-coordinate in first column and y-coordinate in second column
            node = np.reshape(node[:, 0], (len(node), 1))
            node_tag = node_tag - np.ones(len(node_tag))
            node_tag = node_tag.astype('int')
            node[node_tag, :] = node

            # Getting edges (i.e. elements in line_mesh)
            ind_elems = np.where(element_types == 1)[0]  # Indices of elements
            elements = np.array(node_tags_elements[ind_elems[0]])
            num_elements = int(len(elements) / 2)  # Number of elements
            elem_to_node = np.reshape(
                elements, (num_elements, 2)) - np.ones((num_elements, 1))
            elem_to_node = elem_to_node.astype('int')

            node_to_regi = Geometry.GmshMeshReader.get_physical_groups_per_element(
                0, node, elem_to_node)
            elem_to_regi = Geometry.GmshMeshReader.get_physical_groups_per_element(
                1, node, elem_to_node)

            return (node, elem_to_node, node_to_regi, elem_to_regi)

        @staticmethod
        def read_tet_mesh():
            """Returns a tetrahedral pyrit msh model for the given geometry.

            Unifies the generated gmsh model (typically generated via Geometry.create_mesh) with
            the physical groups defined for this Geometry object into one TetMesh object.

            Returns
            -------
            msh : TetMesh
                Instance of a pyrit tetrahedral mesh.

            """
            node_tag, node, _ = gmsh.model.mesh.getNodes()
            element_types, element_tags, node_tags_elements = gmsh.model.mesh.getElements()
            if len(element_types) == 0:
                raise Exception("Cannot read mesh because no elements exist. Have you created the mesh?")

            # Getting tags
            ind_points = np.where(element_types == 15)[0][0]  # Indices of points. 15 is the element type of points
            ind_edges = np.where(element_types == 1)[0][0]  # Indices of edges. 1 is the element type of an edge
            ind_faces = np.where(element_types == 2)[0][0]  # Indices of triangles. 2 is the
            ind_elements = np.where(element_types == 4)[0][0]  # Indices of tetrahedrals

            # element type of an triangle
            point_tags = element_tags[ind_points]
            edge_tags = element_tags[ind_edges]
            triangle_tags = element_tags[ind_faces]
            tetrahedral_tags = element_tags[ind_elements]

            # Getting nodes (points)
            num_node = int(len(node) / 3)  # Number of nodes
            node = np.reshape(node, (num_node, 3))
            node = node[:, 0:3]
            node_tag = node_tag - np.ones(len(node_tag))
            node_tag = node_tag.astype('int')
            node[node_tag, :] = node  # TODO: is this still correct? different in read_trimesh

            # Getting edges (lines)
            edges = np.array(node_tags_elements[ind_edges])
            num_edges = int(len(edges) / 2)  # Number of edges
            gmsh_edge2node = np.reshape(edges, (num_edges, 2)) - 1
            gmsh_edge2node = gmsh_edge2node.astype('int')
            gmsh_edge2node = np.sort(gmsh_edge2node, axis=1)

            # Getting faces (triangles)
            faces = np.array(node_tags_elements[ind_faces])
            num_face = int(len(faces) / 3)  # Number of faces
            gmsh_face2node = np.reshape(faces, (num_face, 3)) - np.ones((num_face, 1))
            gmsh_face2node = gmsh_face2node.astype('int')

            # Getting elements (tetrahedrals)
            elements = np.array(node_tags_elements[ind_elements])
            num_elem = int(len(elements) / 4)  # Number of elements
            elem_to_node = np.reshape(elements, (num_elem, 4)) - np.ones((num_elem, 1))
            # Connection between elements and nodes. Each line contains the indices of the contained nodes
            elem_to_node = elem_to_node.astype('int')  # same for gmsh and pyrit

            msh = TetMesh(node, elem_to_node)
            Geometry.GmshMeshReader.remove_unused(msh, TriMesh.calc_node2elem(elem_to_node, num_node))

            # Physical Groups of nodes
            pg = gmsh.model.getPhysicalGroups(dim=0)
            msh.node2regi = -1 * np.ones(msh.num_node, dtype=int)
            for k, _ in enumerate(pg):
                node_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind = np.array([np.where(point_tags == m)[0] for m in node_tags_pg]).flatten()
                msh.node2regi[ind] = pg[k][1]

            # Physical groups of edges
            pg = gmsh.model.getPhysicalGroups(dim=1)
            msh.edge2regi = -1 * np.ones(msh.num_edge, dtype=int)
            for k, _ in enumerate(pg):
                edge_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind_tmp = np.array([np.where(edge_tags == m)[0] for m in edge_tags_pg]).flatten()
                tmp = gmsh_edge2node[ind_tmp, :]  # nodes for this pg in gmsh indices
                ind = np.array([np.where(np.all(msh.edge2node[:, :] == tmp[m, :], 1))[0][0] for m in
                                range(len(tmp))])  # node indices in pyrit for this pg
                msh.edge2regi[ind] = pg[k][1]

            # Physical groups of faces
            pg = gmsh.model.getPhysicalGroups(dim=2)
            msh.face2regi = -1 * np.ones(msh.num_face, dtype=int)
            for k, _ in enumerate(pg):
                # find tags of faces in gmsh with the physical group
                face_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                # find the indices in the gmsh list of faces corresponding to these tags
                ind_tmp = np.array([np.where(triangle_tags == m)[0] for m in face_tags_pg]).flatten()
                # find the face-to-node entries for these faces. Note the sorting neglects
                # orientation
                tmp = np.sort(gmsh_face2node[ind_tmp, :], 1)
                # transform from gmsh-indices to pyrit by finding equivalent faces
                ind = np.array([np.where(np.all(msh.face2node[:, :] == tmp[m, :], 1))[0][0] for m in range(len(tmp))])
                msh.face2regi[ind] = pg[k][1]

            # Physical groups of tetrahedrals
            pg = gmsh.model.getPhysicalGroups(dim=3)
            msh.elem2regi = -1 * np.ones(num_elem, dtype=int)
            for k, _ in enumerate(pg):
                elem_tags_pg = Geometry.GmshMeshReader.get_pg_per_element(pg[k][0], pg[k][1])
                ind = np.array([np.where(tetrahedral_tags == m)[0] for m in elem_tags_pg]).flatten()
                msh.elem2regi[ind] = pg[k][1]

            return msh

        @staticmethod
        def get_pg_per_element(dim: int, group_tag):
            tags = gmsh.model.getEntitiesForPhysicalGroup(dim, group_tag)

            elem_tags = np.array([], dtype=int)
            for tag in tags:
                _, elem_tags_tmp, _ = gmsh.model.mesh.getElements(dim, tag)
                elem_tags_tmp = elem_tags_tmp[0].astype('int')
                elem_tags = np.unique(np.concatenate((elem_tags, elem_tags_tmp)))

            return elem_tags

        @staticmethod
        def get_physical_groups_per_element(dim: int, node: np.ndarray, elem_to_node: np.ndarray):
            """
            Returns a np.ndarray that assigns to every entity of dimension dim the physical
            group, if possible. It assigns the ID of the physical group. If for an entity
            no physical group exists, -1 is assigned

            Parameters
            ----------
            dim : int
                The dimension of the entity and the physical group.
            node : np.ndarray
                Array containing the coordinates of all nodes .
            elem_to_node : np.ndarray
                Array that maps each element to the nodes that define it. The number of rows
                is the number of elements and the number of columns is the number of nodes
                that define the element (i.e. dim+1)

            Returns
            -------
            np.ndarray
                (N,)-Array with the ID of a physical group for each entity

            """
            pg = gmsh.model.getPhysicalGroups(dim=dim)

            # Physical Groups
            physical_groups = dict()
            physical_group_names = dict()
            for group in pg:
                physical_groups[group[1]], _ = gmsh.model.mesh.getNodesForPhysicalGroup(
                    group[0], group[1])
                physical_groups[group[1]] = physical_groups[group[1]].astype(
                    'int') - np.ones(len(physical_groups[group[1]]), dtype=int)

                physical_group_names[group[1]] = gmsh.model.getPhysicalName(
                    group[0], group[1])

            if dim == 0:  # Physical groups are points
                node_to_regi = -1 * np.ones(len(node), dtype=np.int32)
                for pg_tag in physical_groups.keys():
                    node_to_regi[physical_groups[pg_tag]] = pg_tag
                return node_to_regi
            if dim == 1:
                elem_to_regi = -1 * \
                               np.ones(elem_to_node.shape[0], dtype=np.int32)
                for pg_tag in physical_groups:
                    for k in range(len(elem_to_regi)):
                        if elem_to_node[k, 0] in physical_groups[pg_tag] and \
                                elem_to_node[k, 1] in physical_groups[pg_tag]:
                            elem_to_regi[k] = pg_tag
                return elem_to_regi
            if dim == 2:
                elem_to_regi = -1 * \
                               np.ones(elem_to_node.shape[0], dtype=np.int32)
                for pg_tag in physical_groups:
                    for k in range(len(elem_to_regi)):
                        if elem_to_node[k, 0] in physical_groups[pg_tag] and \
                                elem_to_node[k, 1] in physical_groups[pg_tag] and \
                                elem_to_node[k, 2] in physical_groups[pg_tag]:
                            elem_to_regi[k] = pg_tag
                return elem_to_regi
            if dim == 3:
                elem_to_regi = -1 * \
                               np.ones(elem_to_node.shape[0], dtype=np.int32)
                for pg_tag in physical_groups:
                    for k in range(len(elem_to_regi)):
                        if elem_to_node[k, 0] in physical_groups[pg_tag] and \
                                elem_to_node[k, 1] in physical_groups[pg_tag] and \
                                elem_to_node[k, 2] in physical_groups[pg_tag] and \
                                elem_to_node[k, 3] in physical_groups[pg_tag]:
                            elem_to_regi[k] = pg_tag
                return elem_to_regi

        @njit
        def edges_of_triangle(ind_triangle, edge_to_node, triangle_to_node):
            """
            Gives the indices (with correct sign) of the edges for the triangle with index
            ind_triangle.

            Parameters
            ----------
            ind_triangle : int
                Index of the triangle (i.e. row in triangle_to_node) to get the edges from.
            edge_to_node : np.ndarray
                (number of edges,2)-ndarray with the indices of two nodes defining the edge
                in every row.
            triangle_to_node : np.ndarray
                (number of triangles,3)-ndarray with the indices of three nodes defining
                the triangle in every row

            Raises
            ------
            Warning
                If something goes wrong.

            Returns
            -------
            out : List[int]
                List of the indices of the edges.

            """
            inds = []
            # Get for every point of the triangle all edges that use this point
            for k in range(3):
                inds.append(np.where(edge_to_node ==
                                     triangle_to_node[ind_triangle, k])[0])

            # inds = [np.where(edge_to_node == triangle_to_node[ind_triangle, k])[0] for k in range(3)]
            out = list()  # np.zeros(3, dtype=np.int64)
            for n1, n2 in zip([0, 1, 2], [1, 2, 0]):
                ind = np.intersect1d(inds[n1], inds[n2])
                if triangle_to_node[ind_triangle, n1] == edge_to_node[ind, 0]:
                    out.append(ind)
                elif triangle_to_node[ind_triangle, n1] == edge_to_node[ind, 1]:
                    out.append(-ind)
                else:
                    raise Warning('Not possible')
            return out

        @njit
        def faces_of_tetrahedron(ind_tet, face_to_node, tet_to_node):
            inds = []
            # Get for every point of tetrahedron all faces that use this point
            for k in range(4):
                inds.append(np.where(face_to_node ==
                                     tet_to_node[ind_tet, k])[0])

            out = list()  # np.zeros(3, dtype=np.int64)
            for n1, n2, n3 in zip([0, 0, 0, 3], [1, 1, 2, 1], [2, 3, 3, 2]):
                ind = np.intersect1d(inds[n1], np.intersect1d(inds[n2], inds[n3]))
                out.append(ind)
            return out

        @njit
        def edges_of_tetrahedron(ind_tet, edge_to_node, tet_to_node):
            inds = []
            # Get for every point of tetrahedron all edges that use this point
            for k in range(4):
                inds.append(np.where(edge_to_node ==
                                     tet_to_node[ind_tet, k])[0])

            out = list()  # np.zeros(3, dtype=np.int64)
            for n1, n2 in zip([0, 0, 0, 1, 2, 3], [1, 2, 3, 2, 3, 1]):
                ind = np.intersect1d(inds[n1], inds[n2])
                out.append(ind)
            return out

        @staticmethod
        def remove_unused(mesh: Mesh, node2elem: List[List[int]]) -> NoReturn:
            """Removes these nodes in a mesh that don't belong to any element.

            Parameters
            ----------
            mesh : Mesh
                Instance of a mesh.
            node2elem : List[List[int]]
                A list of list of nodes for each element. See :py:meth:`~pyrit.mesh.Mesh.calc_node2elem`.
            """
            unused_nodes = [k for k, n_to_e in enumerate(node2elem) if len(n_to_e) == 0]
            if len(unused_nodes) == 0:
                return
            used_nodes = np.setdiff1d(np.arange(mesh.num_node), unused_nodes)
            mapping = -1 * np.ones(mesh.num_node, dtype=int)
            mapping[used_nodes] = np.arange(mesh.num_node - len(unused_nodes))

            elem_to_node = mapping[mesh.elem2node]
            mesh.node = mesh.node[used_nodes, :]
            mesh.elem2node = elem_to_node
