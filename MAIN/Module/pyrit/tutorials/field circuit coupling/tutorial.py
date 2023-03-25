# coding=utf-8
"""Field circuit coupling tutorial

To run this tutorial, the package `schemdraw` needs to be installed. It can be installed with

>>> pip install schemdraw
"""

from dataclasses import dataclass
import numpy as np
import gmsh
import scipy.constants
import scipy.sparse.linalg as splinalg
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib

from pyrit.excitation import Excitations, FieldCircuitCoupling, SolidConductor, StrandedConductor, DistributedSolidConductor
from pyrit.geometry import Geometry, Surface, Rectangle
from pyrit.bdrycond import BdryCond, BCDirichlet
from pyrit.material import Materials, Mat, Conductivity, Reluctivity, Resistivity
from pyrit.mesh import AxiMesh
from pyrit.shapefunction import TriAxisymmetricEdgeShapeFunction
from pyrit.region import Regions
from pyrit.problem import MagneticProblemAxiHarmonic

from pyrit.toolbox.geolibrary.Yoke import Yoke


@dataclass
class CoilInYokeStranded:
    radius_yoke: float = 2e-2
    height_yoke: float = 4e-2
    thickness_yoke: float = 5e-3

    width_coil: float = 5e-3  # Width of the coil
    height_coil: float = 2e-2  # Height of the coil
    mid_coil_r: float = 8e-3  # r-coordinate of the middle point of the coil
    mid_coil_z: float = 2e-2  # z-coordinate of the middle point of the coil

    mu_r_iron: float = 10
    fill_factor: float = 1

    mu_0 = scipy.constants.mu_0
    nu_0 = 1 / mu_0

    def create_geometry_stranded(self, stranded_conductor, mesh_size=0.2):
        geo = Geometry('coil in yoke', show_gui=False, mesh_size_factor=mesh_size)

        air = Mat('air', Reluctivity(self.nu_0), Conductivity(0), Resistivity(0))
        copper = Mat('copper', Reluctivity(0.99 * self.nu_0), Conductivity(0), Resistivity(1 / self.fill_factor / 6e7))
        iron = Mat('iron', Reluctivity(1 / self.mu_r_iron * self.nu_0), Conductivity(1e7), Resistivity(1e-7))
        materials = Materials(air, copper, iron)

        bc_homogen = BCDirichlet(0)
        boundary_conditions = BdryCond(bc_homogen)

        excitations = Excitations(stranded_conductor)

        conductor = geo.create_physical_group(1, 2, f'stranded_conductor')
        filling = geo.create_physical_group(10, 2, 'filling')
        yoke_pg = geo.create_physical_group(11, 2, 'yoke')
        bound = geo.create_physical_group(12, 1, 'bound')

        geo.add_material_to_physical_group(air, filling)
        geo.add_material_to_physical_group(copper, conductor)
        geo.add_material_to_physical_group(iron, yoke_pg)
        geo.add_boundary_condition_to_physical_group(bc_homogen, bound)
        geo.add_excitation_to_physical_group(stranded_conductor, conductor)

        with geo.geometry():
            yoke = Yoke('yoke', width=self.radius_yoke, height=self.height_yoke, thickness=self.thickness_yoke,
                        thickness_inner=self.thickness_yoke / 2, air_gap_thickness=0, anchor=Yoke.Anchor.SOUTH_WEST)
            yoke.add_to_gmsh()

            coil = Rectangle(self.mid_coil_r - self.width_coil / 2, self.mid_coil_z - self.height_coil / 2, 0,
                             self.width_coil, self.height_coil)
            coil.add_to_gmsh()

            surf = Surface(yoke.inner_lines, coil.ext_curves)
            surf.add_to_gmsh()

            conductor.add_entity(coil)
            filling.add_entity(surf)
            yoke_pg.add_entity(yoke.surface)
            bound.add_entity(*yoke.outer_lines)

            gmsh.model.mesh.field.add("Box", 7)
            gmsh.model.mesh.field.setNumber(7, "VIn", 0.01 * mesh_size)
            gmsh.model.mesh.field.setNumber(7, "VOut", mesh_size)
            gmsh.model.mesh.field.setNumber(7, "XMin", 0)
            gmsh.model.mesh.field.setNumber(7, "XMax", self.thickness_yoke)
            gmsh.model.mesh.field.setNumber(7, "YMin", 0)
            gmsh.model.mesh.field.setNumber(7, "YMax", self.height_yoke)
            gmsh.model.mesh.field.setAsBackgroundMesh(7)

            geo.create_mesh(2)
            mesh = geo.get_mesh(2)
            regions = geo.get_regions()

        mesh = AxiMesh.generate_axi_mesh(mesh)
        sft = TriAxisymmetricEdgeShapeFunction(mesh)

        return MagneticProblemAxiHarmonic('Coil in Yoke, Stranded conductor', mesh, sft, regions, materials,
                                          boundary_conditions, excitations, frequency=1)


@dataclass
class CoilInYokeSolid:
    radius_yoke: float = 2e-2
    height_yoke: float = 4e-2
    thickness_yoke: float = 5e-3

    width_coil: float = 5e-3  # Width of the coil
    height_coil: float = 2e-2  # Height of the coil
    mid_coil_r: float = 8e-3  # r-coordinate of the middle point of the coil
    mid_coil_z: float = 2e-2  # z-coordinate of the middle point of the coil

    radius_conductor: float = 5e-4  # Radius of one strand
    distance_conductors: float = 1e-4  # Distance from a conductor to its neighbor and to the boundary of the coil
    windings: int = 20

    mu_r_iron: float = 10

    mu_0 = scipy.constants.mu_0
    nu_0 = 1 / mu_0

    @property
    def R_ana(self):
        self.update_distance_conductors()
        rho = 1 / 6e7
        area = self.radius_conductor ** 2 * np.pi
        r_mid, _ = self.compute_mid_points()
        length = sum([2 * r * np.pi for r in r_mid])
        return rho * length / area

    @property
    def L_ana(self):
        area = self.mid_coil_r ** 2 * np.pi
        return self.windings ** 2 * self.mu_0 * area / self.height_coil

    @property
    def fill_factor(self):
        area = self.width_coil * self.height_coil
        area_conductors = self.windings * self.radius_conductor ** 2 * np.pi
        return area_conductors / area

    def check_r_d(self):
        """check if the radius and the distance match."""
        k = (self.height_coil - self.distance_conductors) / (2 * self.radius_conductor + self.distance_conductors)
        if np.isclose(k, np.round(k)):
            return np.round(k)
        return None

    def update_distance_conductors(self):
        k = np.round(
            (self.height_coil - self.distance_conductors) / (2 * self.radius_conductor + self.distance_conductors))
        self.distance_conductors = (self.height_coil - k * 2 * self.radius_conductor) / (k + 1)
        return self.distance_conductors

    @property
    def max_number_windings(self):
        k = (self.height_coil - self.distance_conductors) / (2 * self.radius_conductor + self.distance_conductors)
        m = np.floor(
            (self.width_coil - self.distance_conductors) / (2 * self.radius_conductor + self.distance_conductors))
        return k * m

    def compute_mid_points(self):
        if self.windings > self.max_number_windings:
            raise ValueError("Too many windings")

        num_row = self.check_r_d()

        full_columns = int(np.floor(self.windings / num_row))
        num_in_part_column = int(np.mod(self.windings, num_row))

        coil_bottom_left_r = self.mid_coil_r - self.width_coil / 2
        coil_bottom_left_z = self.mid_coil_z - self.height_coil / 2
        coil_top_left_z = self.mid_coil_z + self.height_coil / 2
        r_coordinates = np.concatenate([np.repeat(
            coil_bottom_left_r + self.radius_conductor + self.distance_conductors + k * (
                    2 * self.radius_conductor + self.distance_conductors), num_row) for k in range(full_columns)])

        if num_in_part_column > 0:
            r_coordinates = np.concatenate([r_coordinates,
                                            np.repeat(
                                                coil_bottom_left_r + self.radius_conductor + self.distance_conductors + full_columns * (
                                                        2 * self.radius_conductor + self.distance_conductors),
                                                num_in_part_column)])

        z_coordinates = np.tile(np.arange(start=coil_bottom_left_z + self.radius_conductor + self.distance_conductors,
                                          stop=coil_top_left_z - 0.99 * self.radius_conductor - self.distance_conductors,
                                          step=2 * self.radius_conductor + self.distance_conductors), full_columns)

        if num_in_part_column > 0:
            b = (self.height_coil - 2 * num_in_part_column * self.radius_conductor - (
                    num_in_part_column - 1) * self.distance_conductors) / 2
            z_coordinates = np.concatenate(
                [z_coordinates, np.arange(start=coil_bottom_left_z + b + self.radius_conductor,
                                          stop=coil_top_left_z - 0.99 * b - self.radius_conductor,
                                          step=2 * self.radius_conductor + self.distance_conductors)])

        return r_coordinates, z_coordinates

    def create_geometry(self, mesh_size=0.2, mesh_size_conductors=0.05):
        self.update_distance_conductors()
        geo = Geometry('coil in yoke', show_gui=False, mesh_size_factor=mesh_size)

        air = Mat('air', Reluctivity(self.nu_0), Conductivity(0), Resistivity(0))
        copper = Mat('copper', Reluctivity(0.99 * self.nu_0), Conductivity(6e7), Resistivity(1 / 6e7))
        iron = Mat('iron', Reluctivity(1 / self.mu_r_iron * self.nu_0), Conductivity(1e7), Resistivity(1e-7))
        insulation = Mat('insulation', Reluctivity(self.nu_0), Conductivity(0), Resistivity(0))
        materials = Materials(air, copper, iron, insulation)

        bc_homogen = BCDirichlet(0)
        boundary_conditions = BdryCond(bc_homogen)
        # dsc = DistributedSolidConductor(self.windings, current=1)
        fccs = [SolidConductor(current=1) for _ in range(self.windings)]
        excitations = Excitations(*fccs)
        # excitations = Excitations(dsc)
        filling = geo.create_physical_group(10, 2, 'filling')
        yoke_pg = geo.create_physical_group(11, 2, 'yoke')
        bound = geo.create_physical_group(12, 1, 'bound')
        insulation_pg = geo.create_physical_group(13, 2, 'insulation')

        conductors = [geo.create_physical_group(100 + k, 2, f'stranded_conductor_{k}') for k in range(self.windings)]

        geo.add_material_to_physical_group(air, filling)
        geo.add_material_to_physical_group(iron, yoke_pg)
        geo.add_material_to_physical_group(insulation, insulation_pg)
        geo.add_boundary_condition_to_physical_group(bc_homogen, bound)
        for k in range(self.windings):
            geo.add_excitation_to_physical_group(fccs[k], conductors[k])
            geo.add_material_to_physical_group(copper, conductors[k])

        with geo.geometry():
            yoke = Yoke('yoke', width=self.radius_yoke, height=self.height_yoke, thickness=self.thickness_yoke,
                        thickness_inner=self.thickness_yoke / 2, air_gap_thickness=0,
                        mesh_size=mesh_size_conductors * 5, anchor=Yoke.Anchor.SOUTH_WEST)
            yoke.add_to_gmsh()

            coil = Rectangle(self.mid_coil_r - self.width_coil / 2, self.mid_coil_z - self.height_coil / 2, 0,
                             self.width_coil, self.height_coil)
            coil.add_to_gmsh()

            surf = Surface(yoke.inner_lines, coil.ext_curves)
            surf.add_to_gmsh()

            conductor_surfaces = []
            r_coordinates, z_coordinates = self.compute_mid_points()
            for k in range(self.windings):
                r = r_coordinates[k]
                z = z_coordinates[k]
                conductor_surface = gmsh.model.occ.addDisk(r, z, 0, self.radius_conductor, self.radius_conductor)
                conductor_surfaces.append(conductor_surface)
                conductors[k].add_entity(conductor_surface)

            # gmsh.model.occ.cut([(2, coil), ], [(2, conductor_surfaces[0]), ], removeTool=False)
            dimtag, _ = gmsh.model.occ.cut([(2, coil.tag), ], [(2, cs) for cs in conductor_surfaces], removeTool=False)

            filling.add_entity(surf)
            yoke_pg.add_entity(yoke.surface)
            bound.add_entity(*yoke.outer_lines)
            insulation_pg.add_entity(dimtag[0][1])

            gmsh.model.mesh.field.add("Box", 7)
            gmsh.model.mesh.field.setNumber(7, "VIn", mesh_size_conductors)
            gmsh.model.mesh.field.setNumber(7, "VOut", mesh_size)
            gmsh.model.mesh.field.setNumber(7, "XMin", self.mid_coil_r - self.width_coil / 2)
            gmsh.model.mesh.field.setNumber(7, "XMax", self.mid_coil_r + self.width_coil / 2)
            gmsh.model.mesh.field.setNumber(7, "YMin", self.mid_coil_z - self.height_coil / 2)
            gmsh.model.mesh.field.setNumber(7, "YMax", self.mid_coil_z + self.height_coil / 2)
            # gmsh.model.mesh.field.setAsBackgroundMesh(7)

            gmsh.model.mesh.field.add("Box", 8)
            gmsh.model.mesh.field.setNumber(8, "VIn", 0.005 * mesh_size)
            gmsh.model.mesh.field.setNumber(8, "VOut", mesh_size)
            gmsh.model.mesh.field.setNumber(8, "XMin", 0)
            gmsh.model.mesh.field.setNumber(8, "XMax", self.thickness_yoke)
            gmsh.model.mesh.field.setNumber(8, "YMin", 0)
            gmsh.model.mesh.field.setNumber(8, "YMax", self.height_yoke)
            # gmsh.model.mesh.field.setAsBackgroundMesh(8)

            gmsh.model.mesh.field.add("Min", 9)
            gmsh.model.mesh.field.setNumbers(9, "FieldsList", [7, 8])

            gmsh.model.mesh.field.setAsBackgroundMesh(9)

            geo.create_mesh(2)
            mesh = geo.get_mesh(2)
            regions = geo.get_regions()

        mesh = AxiMesh.generate_axi_mesh(mesh)
        sft = TriAxisymmetricEdgeShapeFunction(mesh)

        return MagneticProblemAxiHarmonic('Coil in Yoke, Solid coductor', mesh, sft, regions, materials,
                                          boundary_conditions, excitations, frequency=1)


def simulate_stranded():
    coil_in_yoke = CoilInYokeStranded(mu_r_iron=1000, fill_factor=0.565486)
    sc = StrandedConductor(voltage=1e-3, windings=72)
    problem = coil_in_yoke.create_geometry_stranded(sc, mesh_size=0.1)

    problem.frequency = 5e1
    shape_function: TriAxisymmetricEdgeShapeFunction = problem.shape_function
    mesh: AxiMesh = problem.mesh

    curlcurl = shape_function.curlcurl_operator(problem.regions, problem.materials, Reluctivity)
    mass = shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)
    mass_one = shape_function.mass_matrix(1)
    mass_rho = shape_function.mass_matrix(problem.regions, problem.materials, Resistivity)
    rhs = shape_function.load_vector(0)
    matrix = curlcurl + 1j * problem.omega * mass

    matrix_shrink, rhs_shrink = FieldCircuitCoupling.shrink(matrix, rhs, problem, integration_order=2)
    solution_shrink = splinalg.spsolve(matrix_shrink.tocsr(), rhs_shrink.tocsr())
    solution, fcc_ids = FieldCircuitCoupling.inflate(solution_shrink, problem)

    current = sc.current
    voltage = sc.voltage
    impedance = sc.impedance
    resistance = sc.resistance
    inductance = sc.inductance(problem.omega)

    print(f"Resistance: {resistance}\tInductance: {inductance}")

    time = np.linspace(0, 3 / problem.frequency, 200)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(time, sc.voltage * np.cos(problem.omega * time), label='Voltage')
    ax.plot(time, np.abs(current) * np.cos(problem.omega * time + np.angle(current)), label="current")
    ax.legend()
    plt.show()

    nodes_not_on_axis = np.setdiff1d(np.arange(mesh.num_node), mesh.nodes_on_axis)
    # mesh.plot_equilines(np.real(solution / (2 * np.pi)), 20)

    b = shape_function.curl(solution)
    mesh.plot_scalar_field(np.hypot(np.real(b[:, 0]), np.real(b[:, 1])), title="flux density")

    e_eddy = -1j * problem.omega * solution
    j_source = None
    for fcc in StrandedConductor.iter_fcc(problem.excitations):
        fcc.compute_current_distribution_function(problem)
        if j_source is None:
            j_source = fcc.current * fcc.windings * fcc.current_distribution_function
        else:
            j_source = j_source + fcc.current * fcc.windings * fcc.current_distribution_function

    e_source = project_to_c_vector(mass_rho @ j_source, mass_one.tocsr(), mesh.num_node, nodes_not_on_axis)

    e = e_eddy + e_source
    j = j_source.toarray().flatten() + project_to_c_vector(mass @ e, mass_one.tocsr(), mesh.num_node, nodes_not_on_axis)
    e_field = field_vector(e, mesh, nodes_not_on_axis)
    j_field = field_vector(j, mesh, nodes_not_on_axis)

    mesh.plot_scalar_field(np.real(e_field), title="e field")
    mesh.plot_scalar_field(np.real(j_field), title="j field")
    plt.show()


def field_vector(vector, mesh, nodes_not_on_axis):
    field = np.zeros_like(vector, dtype=vector.dtype)
    field[nodes_not_on_axis] = 1 / (2 * np.pi * mesh.node[nodes_not_on_axis, 0]) * vector[nodes_not_on_axis]
    return field


def project_to_c_vector(r_vector, matrix_one, num_node, nodes_not_on_axis):
    matrix_shrink = matrix_one[np.ix_(nodes_not_on_axis, nodes_not_on_axis)]
    s_vector_shrink = r_vector[nodes_not_on_axis]

    tmp = splinalg.spsolve(matrix_shrink, s_vector_shrink)
    c_vector = np.zeros(num_node, dtype=r_vector.dtype)
    c_vector[nodes_not_on_axis] = tmp
    return c_vector


def simulate_solid():
    coil_in_yoke = CoilInYokeSolid(windings=72, mu_r_iron=1000)
    problem = coil_in_yoke.create_geometry(mesh_size=0.2, mesh_size_conductors=0.0005)

    print(f"Fill_factor: {coil_in_yoke.fill_factor}")
    print(f"Analytic R: {coil_in_yoke.R_ana}")
    print(f"Analytic L: {coil_in_yoke.L_ana}")
    problem.frequency = 5e4
    shape_function: TriAxisymmetricEdgeShapeFunction = problem.shape_function
    mesh: AxiMesh = problem.mesh

    curlcurl = shape_function.curlcurl_operator(problem.regions, problem.materials, Reluctivity)
    mass = shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)
    mass_one = shape_function.mass_matrix(1)
    rhs = shape_function.load_vector(0)
    matrix = curlcurl + 1j * problem.omega * mass

    matrix_shrink, rhs_shrink = FieldCircuitCoupling.shrink(matrix, rhs, problem, integration_order=2)
    solution_shrink = splinalg.spsolve(matrix_shrink.tocsr(), rhs_shrink.tocsr())
    solution, ids = FieldCircuitCoupling.inflate(solution_shrink, problem)

    voltage = sum((problem.excitations.get_exci(id).voltage for id in ids))
    print(voltage)

    time = np.linspace(0, 3 / problem.frequency, 200)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(time, np.cos(problem.omega * time), label='Current')
    ax.plot(time, np.abs(voltage) * np.cos(problem.omega * time + np.angle(voltage)), label="voltage")
    ax.legend()
    # plt.show()

    # Resistance and Inductance
    impedance = voltage / 1
    resistance = impedance.real
    inductance = impedance.imag / problem.omega

    print(f"Resistance: {resistance}\tInductance: {inductance}")

    nodes_not_on_axis = np.setdiff1d(np.arange(mesh.num_node), mesh.nodes_on_axis)
    mesh.plot_equilines(np.real(solution / (2 * np.pi)), 20, colors='k')

    b = shape_function.curl(solution)
    mesh.plot_scalar_field(np.hypot(np.real(b[:, 0]), np.real(b[:, 1])), title="Magnetic flux density")

    e_eddy = -1j * problem.omega * solution

    e_source = sparse.coo_matrix((mesh.num_node, 1))
    for key in ids:
        # noinspection PyUnresolvedReferences
        fcc: SolidConductor = problem.excitations.get_exci(key)
        fcc.compute_voltage_distribution_function(problem)
        e_source = e_source + fcc.voltage * fcc.voltage_distribution_function

    e = e_eddy + e_source.toarray().flatten()
    j = project_to_c_vector(mass @ e, mass_one.tocsr(), mesh.num_node, nodes_not_on_axis)
    e_field = field_vector(e, mesh, nodes_not_on_axis)
    j_field = field_vector(j, mesh, nodes_not_on_axis)
    mesh.plot_scalar_field(np.real(e_field), title="Electric field")
    mesh.plot_scalar_field(np.real(j_field), title="Current density")
    mesh.plot_scalar_field(np.real(solution), title="vector potential")
    plt.show()


# noinspection PyPackageRequirements
def draw_circuit():
    import schemdraw
    import schemdraw.elements as elm

    with schemdraw.Drawing() as d:
        d += elm.SourceSin().up().label(r'$V_s$=1V')
        d += (R := elm.Resistor().right().label('R'))
        d += elm.CurrentLabelInline().at(R).label(r'$I_s$')
        d += elm.CurrentLabel(top=False).at(R).label(r'$V_R$')
        d += elm.Dot()
        d.push()
        d += elm.Line().right()
        d += (L := elm.Inductor().down())
        d += elm.CurrentLabelInline().at(L).label(r'$I_L$')
        d += elm.CurrentLabel(top=False, reverse=True).at(L).label('V')
        d += elm.Line().left()
        d += elm.Dot()
        d += elm.Line().left()
        d.pop()
        d += (C := elm.Capacitor().down().label('C'))
        d += elm.CurrentLabelInline().at(C).label(r'$I_C$', loc='bottom')


def simulate_circuit_coupling_solid():
    """Couple the coil to an external circuit with a resistance and a capacitor

    Returns
    -------

    """
    resistance = 2  # in ohm
    capacitance = 1e-6  # in farad

    windings = 72
    coil_in_yoke = CoilInYokeSolid(windings=windings, mu_r_iron=1000)
    problem = coil_in_yoke.create_geometry(mesh_size=0.2, mesh_size_conductors=0.0009)

    problem.frequency = 5e4
    shape_function: TriAxisymmetricEdgeShapeFunction = problem.shape_function
    mesh: AxiMesh = problem.mesh

    curlcurl = shape_function.curlcurl_operator(problem.regions, problem.materials, Reluctivity)
    mass = shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)
    mass_one = shape_function.mass_matrix(1)
    rhs = shape_function.load_vector(0)
    matrix = curlcurl + 1j * problem.omega * mass

    solid_conductors = SolidConductor.iter_fcc(problem.excitations)

    # region From FieldCircuitCoupling.shrink

    bottom_matrices, right_matrices, diagonal_matrices, rhs_matrices = [], [], [], []

    for fcc in solid_conductors:
        # get the coupling values
        # noinspection PyTupleAssignmentBalance
        fcc.compute_coupling_values(problem)

        bottom_matrices.append(fcc.bottom_coupling_matrix)
        right_matrices.append(fcc.right_coupling_matrix)
        diagonal_matrices.append(sparse.coo_matrix(np.array([[fcc.coupling_value]])))

    # Assembly of the matrices
    bottom_matrix: sparse.csr_matrix = sparse.vstack(bottom_matrices, format='csr')
    right_matrix: sparse.csr_matrix = sparse.hstack(right_matrices, format='csr')
    diagonal_matrix = sparse.block_diag(diagonal_matrices, format='csr')
    rhs_fcc = sparse.coo_matrix((windings, 1))

    # shrink the system with the sft
    ret = problem.shape_function.shrink(matrix, rhs, problem, integration_order=2)
    matrix_shrink_sf, rhs_shrink_sf, indices_not_dof, _, _ = ret

    # compose the coupling matrices and coupling values with the information from the sft.shrink
    all_indices = np.arange(rhs.shape[0])
    # num_fcc = len(bottom_matrices)
    new_lines = bottom_matrix.shape[0]
    new_shape = matrix_shrink_sf.shape[0]
    indices_dof = np.setdiff1d(all_indices, indices_not_dof)
    bottom_matrix = bottom_matrix[np.ix_(np.arange(new_lines), indices_dof)]
    right_matrix = right_matrix[np.ix_(indices_dof, np.arange(new_lines))]
    bottom_matrix.resize((bottom_matrix.shape[0], new_shape))
    right_matrix.resize((new_shape, right_matrix.shape[1]))

    # compute the final system
    matrix_shrink = sparse.bmat([[matrix_shrink_sf, right_matrix],
                                 [bottom_matrix, diagonal_matrix]])
    rhs_shrink = sparse.vstack([rhs_shrink_sf, rhs_fcc])

    # endregion

    # Defining new matrices
    old_size = matrix_shrink.shape[0]
    bottom_matrix = sparse.coo_matrix(
        (-1 * np.ones(windings), (np.zeros(windings), np.arange(old_size - windings, old_size))), shape=(3, old_size))
    right_matrix = bottom_matrix.transpose()
    bottom_right_matrix = np.array([[0, 1, 0],
                                    [1, 1j * problem.omega * capacitance, -1 / resistance],
                                    [0, 1, 1]])
    bottom_right_matrix = sparse.coo_matrix(bottom_right_matrix)

    new_matrix = sparse.bmat([[matrix_shrink, right_matrix],
                              [bottom_matrix, bottom_right_matrix]], format='csr')

    new_rhs = rhs_shrink.tocsr()
    new_rhs = sparse.vstack([new_rhs, sparse.coo_matrix(([1], ([2], [0])))], format='csr')

    solution = splinalg.spsolve(new_matrix, new_rhs)
    voltages = solution[-3-windings:-3]
    I_L = solution[-3]
    V = solution[-2]
    V_R = solution[-1]

    solution = problem.shape_function.inflate(solution[:-3-windings], problem)
    for k, fcc in enumerate(solid_conductors):
        fcc.voltage = voltages[k]
        fcc.activate_voltage_mode()
        fcc.current = I_L

    # Postprocessing
    I_s = V_R / resistance
    I_c = I_s - I_L

    impedance = V / I_L
    resistance = impedance.real
    inductance = impedance.imag / problem.omega
    print(f"Resistance: {resistance}\tInductance: {inductance}")

    print(f"Voltages:\nTotal voltage: {1}\nVoltage over R: {V_R}\nVoltage over L,C: {V}\n")
    print(f"Currents:\nTotal current: {I_s}\nCurrent through L: {I_L}\nCurrent through C: {I_c}\n")
    time = np.linspace(0, 3 / problem.frequency, 200)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(time, np.cos(problem.omega * time), label=r'$V_s$')
    ax.plot(time, np.abs(V_R) * np.cos(problem.omega * time + np.angle(V_R)), label=r"$V_R$")
    ax.plot(time, np.abs(V) * np.cos(problem.omega * time + np.angle(V)), label=r"$V$")
    ax.set_title('Voltages')
    ax.legend()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(time, np.abs(I_s) * np.cos(problem.omega * time + np.angle(I_s)), label=r"$I_s$")
    ax.plot(time, np.abs(I_c) * np.cos(problem.omega * time + np.angle(I_c)), label=r"$I_c$")
    ax.plot(time, np.abs(I_L) * np.cos(problem.omega * time + np.angle(I_L)), label=r"$I_L$")
    ax.set_title('Currents')
    ax.legend()

    plt.show()

    nodes_not_on_axis = np.setdiff1d(np.arange(mesh.num_node), mesh.nodes_on_axis)
    mesh.plot_equilines(np.real(solution / (2 * np.pi)), 20, colors='black')

    b = shape_function.curl(solution)
    mesh.plot_scalar_field(np.hypot(np.real(b[:, 0]), np.real(b[:, 1])), title="Magnetic flux density")

    e_eddy = -1j * problem.omega * solution

    e_source = sparse.coo_matrix((mesh.num_node, 1))
    for fcc in solid_conductors:
        # noinspection PyUnresolvedReferences
        fcc.compute_voltage_distribution_function(problem)
        e_source = e_source + fcc.voltage * fcc.voltage_distribution_function

    e = e_eddy + e_source.toarray().flatten()
    j = project_to_c_vector(mass @ e, mass_one.tocsr(), mesh.num_node, nodes_not_on_axis)
    e_field = field_vector(e, mesh, nodes_not_on_axis)
    j_field = field_vector(j, mesh, nodes_not_on_axis)
    mesh.plot_scalar_field(np.real(e_field), title="Electric field")
    mesh.plot_scalar_field(np.real(j_field), title="Current density")

    plt.show()


def main():
    simulate_solid()

    simulate_stranded()

    draw_circuit()
    simulate_circuit_coupling_solid()


if __name__ == '__main__':
    main()
