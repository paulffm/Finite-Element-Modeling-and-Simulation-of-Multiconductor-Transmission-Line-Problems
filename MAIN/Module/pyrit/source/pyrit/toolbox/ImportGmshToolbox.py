# coding=utf-8
"""
Creates TriMesh from .msh file.

.. sectionauthor:: Ruppert, Bundschuh
"""
# pylint: disable=unspecified-encoding

from typing import Type, Union, TYPE_CHECKING, Dict, Tuple
from pathlib import Path
import numpy as np
import gmsh

from pyrit.excitation import Excitations
from pyrit.geometry import Geometry
from pyrit.material import Materials
from pyrit.bdrycond import BdryCond
from pyrit.region import Regi, Regions

if TYPE_CHECKING:
    from pyrit.mesh import Mesh


def geo_to_msh(path: Union[str, 'Path'], refinement_steps: int = 0, **kwargs) -> Path:
    """Generates .msh file from .geo file.

    Parameters
    ----------
    path: Union[str, 'Path']
        Path of the .geo file (with or without file ending)
    refinement_steps: int, optional
        Number of mesh refinement steps. Default is 0
    kwargs :
        Additional keyword arguments. With 'out_path' the path of the output file can be specified (if different from
        the input file). The remaining keyword arguments are passed to :py:class:`~pyrit.geometry.Geometry.Geometry`.

    Returns
    -------
    out_path : Path
        The path of the .msh file.
    """
    path = Path(path).with_suffix('.geo')
    out_path = kwargs.pop('out_path', path)

    geo = Geometry(path.stem, **kwargs)

    with geo.geometry():
        gmsh.merge(str(path))  # Loading the .geo file

        gmsh.model.mesh.generate(gmsh.model.getDimension())  # Generate mesh

        # Refine the mesh
        for _ in range(refinement_steps):
            geo.refine_mesh()

        gmsh.write(str(out_path.with_suffix('.msh')))  # Write .msh-file

    return out_path


def read_msh_file(path: Union[str, 'Path'], mesh_type: Type['Mesh']) -> 'Mesh':
    """Processes .msh file creates Mesh object.

    Parameters
    ----------
    path: Union[str, 'Path']
        Path of the .msh file (with or without file ending).
    mesh_type : Type[Mesh]
        The class of the returned mesh object.

    Returns
    -------
    mesh: Mesh
        The mesh object.
    """
    path = Path(path).with_suffix('.msh')

    with open(path, 'rt') as f:
        # Read node information
        for line in f:
            # Read node information
            if line == "$Nodes\n":
                num_nodes = int(str.split(f.readline())[0])
                data = np.loadtxt(f, max_rows=num_nodes)
                node = data[:, 1:3]
            if line == "$Elements\n":
                numedgeelem = int(str.split(f.readline())[0])  # number of elements and edges combined
                break

        # Read edge information
        edges_complete = False  # True if all lines containing edge information have been processed
        edge2node = np.empty([0, 2])
        edge2phys = np.empty([0, 1])
        elem2node = np.empty([0, 3])
        elem2phys = np.empty([0, 1])
        count_lines = 0
        while not edges_complete:
            # Read line and increase counter of processed lines
            line = f.readline()
            count_lines += 1

            # Collect space delimited entries of line in list
            data_line = list(map(int, line.split()))

            # If the line has 8 entries, all lines belonging to edges have been processed.
            # This and all following lines then correspond to elements
            if len(data_line) == 8:
                edges_complete = True
                elem2node = np.vstack((elem2node, data_line[5:8]))
                elem2phys = np.vstack((elem2phys, data_line[3]))
            else:
                edge2node = np.vstack((edge2node, data_line[5:7]))
                edge2phys = np.vstack((edge2phys, data_line[3]))

        # Read element data
        number_of_remaining_lines = numedgeelem - count_lines
        data = np.loadtxt(f, max_rows=number_of_remaining_lines)
        elem2node = np.vstack((elem2node, data[:, 5:8]))
        elem2phys = np.squeeze(np.vstack((elem2phys, data[:, 3:4])))

        # index correction because python starts indexing at 0
        elem2node = elem2node - 1
        edge2node = edge2node - 1
        # Convert entries from float to int, s.t. they can be used for indexing later
        elem2node = elem2node.astype(int)
        edge2node = edge2node.astype(int)

    # Create mesh object
    mesh = mesh_type(node, elem2node)

    # mesh.edge2node = edge2node
    mesh.edge2regi = - np.ones((mesh.num_edge,))
    for i in range(0, np.shape(edge2node)[0]):
        idx_edge = np.where(((mesh.edge2node[:, 0] == edge2node[i, 0]) & (mesh.edge2node[:, 1] == edge2node[i, 1]))
                            | ((mesh.edge2node[:, 1] == edge2node[i, 0]) &
                               (mesh.edge2node[:, 0] == edge2node[i, 1])))[0][0]
        mesh.edge2regi[idx_edge] = edge2phys[i, 0].astype(int)

    mesh.elem2regi = elem2phys.astype(int)

    return mesh


def read_physical_groups(path: Union[str, 'Path']) -> Dict[str, Tuple[int]]:
    """Reads information about physical groups from .msh file

    Parameters
    ----------
    path: Union[str, 'Path']
        Path of the .msh file (with or without ending)

    Returns
    -------
    physical_groups : Dict[str, Tuple[int]]
        The key is the name of a physical group as defined in the .msh file. The value is a Tuple of int. The first
        entry is the dimension of the physical group and the second one is its tag.
    """
    path = Path(path).with_suffix('.msh')

    with open(path, 'rt') as f:
        # Read physical group information
        for line in f:
            # Read physical group info
            if line == "$PhysicalNames\n":
                num_pg = int(str.split(f.readline())[0])
                physical_groups = {}
                for _ in range(num_pg):
                    pg_tmp = str.split(f.readline())
                    physical_groups[pg_tmp[2][1:-1]] = (pg_tmp[0], pg_tmp[1])
            if line == "$Nodes\n":
                break

    return physical_groups


def generate_regions(path: Union[str, 'Path'], materials: Materials, bdrycond: BdryCond = None,
                     excitations: Excitations = None) -> Regions:
    """Creates Regions object according to names and tag specified in .msh file.

    The names of the materials and boundary conditions must match the physical groups of the .msh file.

    Parameters
    ----------
    path: Union[str, 'Path']
        Path of the .msh file (with or without file ending)
    materials: Materials
    bdrycond: BdryCond
    excitations : Excitations

    Returns
    -------
    regions: Regions
    """
    physical_groups = read_physical_groups(path)
    regions = Regions()
    for pg_name, pg in physical_groups.items():
        pg_dim = int(pg[0])
        pg_tag = int(pg[1])

        # Check if a material belongs to region
        mat_tag = None
        mat_ids = materials.get_ids()
        mat_names = [materials.get_material(id).name for id in mat_ids]
        if pg_name in mat_names:
            mat_tag = materials.get_material_id_from_name(pg_name)

        # Check if a BC belongs to region
        bc_tag = None
        if bdrycond is not None:
            for bc_id in bdrycond.get_ids():
                bc = bdrycond.get_bc(bc_id)
                if bc.name == pg_name:
                    bc_tag = bc_id

        exci_tag = None
        if excitations is not None:
            for exci_id in excitations.get_ids():
                exci = excitations.get_exci(exci_id)
                if exci.name == pg_name:
                    exci_tag = exci_id
        regions.add_regi(Regi(pg_tag, pg_dim, mat=mat_tag, bc=bc_tag, exci=exci_tag))

    return regions
