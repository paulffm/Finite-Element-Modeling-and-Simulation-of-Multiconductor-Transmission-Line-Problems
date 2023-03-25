import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt
from dataclasses import dataclass
from pyrit.geometry import Geometry, Circle, Surface
from scipy.constants import mu_0, epsilon_0
from pyrit.bdrycond import BCDirichlet, BdryCond, BCFloating
from pyrit.material import Materials, Mat, Permittivity, Conductivity, Reluctivity
from pyrit.shapefunction import TriCartesianNodalShapeFunction
from pyrit.excitation import Excitations, CurrentDensity
from pyrit.problem import ElectricProblemCartStatic
from pyrit.problem import Problem
from pyrit.mesh import TriMesh

show_plot = True


def z2polar(z):
    return (abs(z), np.angle(z, deg=True) )


@dataclass
class PowerCable:
    # given
    wire_radius: float = 1.1e-3
    radius_wires_center_points: float = 1.5e-3
    outer_conductor_inner_radius: float = 3e-3
    outer_conductor_outer_radius: float = 3.2e-3
    conductivity_copper: float = 57.7e6
    conductivity_surrounding: float = 0
    relative_permittivity_surrounding: float = 10

    # Strom?
    current: float = 16
    model_name: str = "PowerCable_ys"
    depth: float = 1

    def __post_init__(self):
        self.id_wire_u = 1
        self.id_wire_v = 2
        self.id_wire_w = 3
        self.id_insulation = 10
        self.id_outer_conductor = 11
        self.id_outer_bound = 12

    @property
    def ids_wires(self):
        return [self.id_wire_u, self.id_wire_v, self.id_wire_w]

    @property
    def current_density(self):
        """Current density in [A/m^2]"""
        return self.current / (np.pi * self.wire_radius ** 2)

    def create_problem(self, bcs, magn=True, **kwargs):
        geo = Geometry("Power cable", **kwargs)

        # fragen: Reluktanz?
        materials = Materials(
            mat_wire_u := Mat("Wire u", Permittivity(epsilon_0), Reluctivity(1 / mu_0),
                              Conductivity(self.conductivity_copper)),
            mat_wire_v := Mat("Wire v", Permittivity(epsilon_0), Reluctivity(1 / mu_0),
                              Conductivity(self.conductivity_copper)),
            mat_wire_w := Mat("Wire w", Permittivity(epsilon_0), Reluctivity(1 / mu_0),
                              Conductivity(self.conductivity_copper)),
            mat_outer_cond := Mat("Outer conductor", Permittivity(epsilon_0), Reluctivity(1 / mu_0),
                                  Conductivity(self.conductivity_copper)),
            mat_insulation := Mat("Insulation", Permittivity(self.relative_permittivity_surrounding * epsilon_0),
                                  Reluctivity(1 / mu_0),
                                  Conductivity(self.conductivity_surrounding)))


        # given
        pg_wire_u = geo.create_physical_group(self.id_wire_u, 2, "Wire u")
        pg_wire_v = geo.create_physical_group(self.id_wire_v, 2, "Wire v")
        pg_wire_w = geo.create_physical_group(self.id_wire_w, 2, "Wire w")
        pg_insulation = geo.create_physical_group(self.id_insulation, 2, "Insulation")
        pg_outer_conductor = geo.create_physical_group(self.id_outer_conductor, 2, "Outer conductor")
        pg_outer_bound = geo.create_physical_group(self.id_outer_bound, 1, "Outer bound")

        max_tag = 1
        number_wires = 3
        # Physical groups of the boundaries of the right conductors
        pgs_conductor_bound = [geo.create_physical_group(max_tag + 12 + k, 1, f"conductor_bound{k}") for k
                               in range(number_wires)]

        # add material
        geo.add_material_to_physical_group(mat_insulation, pg_insulation)
        geo.add_material_to_physical_group(mat_wire_u, pg_wire_u)
        geo.add_material_to_physical_group(mat_wire_v, pg_wire_v)
        geo.add_material_to_physical_group(mat_wire_w, pg_wire_w)
        geo.add_material_to_physical_group(mat_outer_cond, pg_outer_conductor)


        # creating boundary condition
        boundary_conditions = BdryCond(bc_outer := BCDirichlet(0), *bcs)
        for bc, pg in zip(bcs, pgs_conductor_bound):
            geo.add_boundary_condition_to_physical_group(bc, pg)

        geo.add_boundary_condition_to_physical_group(bc_outer, pg_outer_conductor)
        # geo.add_boundary_condition_to_physical_group(bc_outer, pg_outer_bound)



        # given
        # %% Building model in gmsh and creating the mesh
        with geo:
            # creating the wires?
            circle_u = Circle(self.radius_wires_center_points, 0, 0, self.wire_radius)
            circle_v = Circle(-0.5 * self.radius_wires_center_points, np.sqrt(3) / 2 * self.radius_wires_center_points,
                              0, self.wire_radius)
            circle_w = Circle(-0.5 * self.radius_wires_center_points,
                              -1 * np.sqrt(3) / 2 * self.radius_wires_center_points, 0, self.wire_radius)

            circle_outer_conductor_inner = Circle(0, 0, 0, self.outer_conductor_inner_radius)
            circle_outer_conductor_outer = Circle(0, 0, 0, self.outer_conductor_outer_radius)

            outer_conductor = Surface([circle_outer_conductor_outer], [circle_outer_conductor_inner])

            # insulation ist alles bis auf outer conductor inner und die wires
            insulation = Surface([circle_outer_conductor_inner], [circle_u, ], [circle_v, ], [circle_w, ])
            surface_u = Surface([circle_u])
            surface_v = Surface([circle_v])
            surface_w = Surface([circle_w])

            for surf in (outer_conductor, insulation, surface_u, surface_v, surface_w):
                surf.add_to_gmsh()

            pg_wire_u.add_entity(surface_u)
            pg_wire_v.add_entity(surface_v)
            pg_wire_w.add_entity(surface_w)
            pg_insulation.add_entity(insulation)
            pg_outer_conductor.add_entity(outer_conductor)
            pg_outer_bound.add_entity(circle_outer_conductor_outer)

            pgs_conductor_bound[0].add_entity(circle_u)
            pgs_conductor_bound[1].add_entity(circle_v)
            pgs_conductor_bound[2].add_entity(circle_w)

            geo.create_mesh(2)
            mesh = geo.get_mesh(2)
            regions = geo.get_regions()

        # Defining the shape function for solving the FE problem

        prb = Problem("Power Cable", mesh, None, regions, materials, boundary_conditions, excitations=[])

        return prb



def main():
    phi_elec = []
    K_list = []

    n = 3

    for i in range(n):
        voltages = np.zeros(n)
        voltages[i] = 1
        print(voltages)
        bcs = [BCDirichlet(val) if val != 0 else BCFloating() for val in voltages]
        power_cable = PowerCable()
        problem = power_cable.create_problem(mesh_size_factor=0.2, bcs=bcs, magn=False, show_gui=False)
        mesh: TriMesh = problem.mesh
        shape_function = TriCartesianNodalShapeFunction(mesh)
        problem.shape_function = shape_function

        # ValueError: Matrix A is singular, because it contains empty row(s)
        divgrad_e = shape_function.divgrad_operator(problem.regions, problem.materials, Permittivity)
        #mass = shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)

        matrix = divgrad_e

        # habe keine load gegeben: wird standardmäßig auf 0 gesetzt: lösbar durch RB
        load = shape_function.load_vector(0)

        matrix_shrink, rhs_shrink, _, _, support_data = shape_function.shrink(matrix, load, problem, 1)
        a_shrink, _ = type(problem).solve_linear_system(matrix_shrink.tocsr(), rhs_shrink.tocsr())
        phi_potential = shape_function.inflate(a_shrink, problem, support_data)
        phi_real = np.real(phi_potential)

        K_list.append(np.asarray(matrix.toarray()))
        phi_elec.append(np.asarray(phi_potential).reshape(-1, 1))
        #print('pot', solution.potential.shape)

        mesh.plot_scalar_field(phi_real, title='Phi distribution')
        mesh.plot_equilines(phi_real, title='Phi Äquipotentiallinien')
        #mesh.plot_scalar_field(np.linalg.norm(np.linalg.norm(phi_potential)), title="Absolute d field")
        plt.show()


    Ke_arr = np.hstack((K_list[0], K_list[1], K_list[2]))
    phi_arr = np.hstack((phi_elec[0], phi_elec[1], phi_elec[2]))

    print('C', phi_arr.T @ K_list[0] @ phi_arr)



if __name__ == '__main__':
    main()