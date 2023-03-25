import numpy as np
from dataclasses import dataclass
import gmsh
import sys
'''
Vorgegebenes Python file für die Erstellung des meshs
'''
# This implementation is well suited for small meshes.
@dataclass
class Mesh:
    node_tag_data : np.ndarray
    node_data : np.ndarray
    elementTypes : np.ndarray
    element_tags : np.ndarray
    nodeTags_elements : np.ndarray

    # treat nodes
    @property
    def num_node(self):
        # number of nodes
        return int(len(self.node_data ) /3)

    @property
    def node(self):
        # nodes
        node = np.reshape(self.node_data, (self.num_node, 3))
        # coordinates of nodes. x-coordinate in first column
        # and y-coordinate in second column
        node = node[:, 0:2]
        return node

    @property
    def node_tag(self):
        # ID of nodes
        node_tag = self.node_tag_data - np.ones(len(self.node_tag_data))
        node_tag = node_tag.astype('int')
        np.put_along_axis(self.node ,np.c_[node_tag ,node_tag] ,self.node ,axis=0)
        return node_tag

        # treat elements
    @property
    def ind_elements(self):
        # index of elements
        return np.where(self.elementTypes == 2)[0]

    @property
    def elements(self):
        # elements
        return np.array(self.nodeTags_elements[self.ind_elements[0]])

    @property
    def num_elements(self):
        # number of elements
        return int(len(self.elements ) /3)

    @property
    def elem_to_node(self):
        # Associate elements (triangles) and their respective nodes.
        # Connection between elements and nodes.
        # Each line contains the indices of the contained nodes
        elem_to_node = np.reshape(self.elements, (self.num_elements, 3)) - np.ones(
            (self.num_elements, 1))
        elem_to_node = elem_to_node.astype('int')
        return elem_to_node

    @property
    def all_edges(self):
        # Contains all edges
        return np.r_[self.elem_to_node[:, [0, 1]],
                     self.elem_to_node[:, [1, 2]],
                     self.elem_to_node[:, [0, 2]]]

    @property
    def edge_to_node(self):
        # Associate edges and their respective nodes.
        # [1 0 ], [1 0], [2 3] unique:[1 0 ] [2 3]
        return np.unique(np.sort(self.all_edges), axis=0)

    # Import Gmsh File
    # Alternatively, the model can be contructed directly with the python interface directly.

    @staticmethod
    def create():
        model_name = 'wire'
        # Initialize gmsh
        gmsh.initialize(sys.argv[0])
        # gmsh.open(model_name+".msh")
        gmsh.model.add(model_name)
        gmsh.merge(model_name + ".geo")
        gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.8)  # control grid size here.
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)  # MATLAB compatible mesh file format
        gmsh.model.occ.synchronize()
        # gmsh.fltk.run()                                  # uncomment, if you want to inspect the geometry
        gmsh.model.mesh.generate(dim=2)  # 2D mesh
        # gmsh.fltk.run()                                  # uncomment, if you want to inspect the mesh
        # gmsh.write(model_name+".msh")                    # writes the mesh to a file

        # Read mesh data
        node_tag, node, _ = gmsh.model.mesh.getNodes()
        elementTypes, element_tags, nodeTags_elements = gmsh.model.mesh.getElements()

        # Process mesh data
        msh = Mesh(node_tag, node, elementTypes, element_tags, nodeTags_elements)
        return msh