import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt
from dataclasses import dataclass
from pyrit.geometry import Geometry, Circle, Surface
from scipy.constants import mu_0, epsilon_0
from pyrit.bdrycond import BCDirichlet, BdryCond
from pyrit.material import Materials, Mat, Permittivity, Conductivity, Reluctivity
from pyrit.shapefunction import TriCartesianEdgeShapeFunction
from pyrit.excitation import Excitations, CurrentDensity
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

    def create_problem(self, magn=True, **kwargs):
        geo = Geometry("Power cable", **kwargs)

        # fragen: Reluktanz?
        materials = Materials(
            mat_wire_u := Mat("Wire u", Permittivity(epsilon_0),  Reluctivity(1 / mu_0),
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

        # creating boundary condition
        boundary_conditions = BdryCond(bc := BCDirichlet(0))

        #excitations = Excitations(exci := CurrentDensity(self.current_density))

        # given
        pg_wire_u = geo.create_physical_group(self.id_wire_u, 2, "Wire u")
        pg_wire_v = geo.create_physical_group(self.id_wire_v, 2, "Wire v")
        pg_wire_w = geo.create_physical_group(self.id_wire_w, 2, "Wire w")
        pg_insulation = geo.create_physical_group(self.id_insulation, 2, "Insulation")
        pg_outer_conductor = geo.create_physical_group(self.id_outer_conductor, 2, "Outer conductor")
        pg_outer_bound = geo.create_physical_group(self.id_outer_bound, 1, "Outer bound")


        # add material
        geo.add_material_to_physical_group(mat_insulation, pg_insulation)
        geo.add_material_to_physical_group(mat_wire_u, pg_wire_u)
        geo.add_material_to_physical_group(mat_wire_v, pg_wire_v)
        geo.add_material_to_physical_group(mat_wire_w, pg_wire_w)
        geo.add_material_to_physical_group(mat_outer_cond, pg_outer_conductor)

        # add excitation
        if magn:
            vals = (self.current_density).tolist()
            print('Current density:', vals)
            excis = [CurrentDensity(val) for val in vals]
            excitations = Excitations(*(excis))
            geo.add_excitation_to_physical_group(excis[0], pg_wire_u)
            geo.add_excitation_to_physical_group(excis[1], pg_wire_v)
            geo.add_excitation_to_physical_group(excis[2], pg_wire_w)


        # set bc
        geo.add_boundary_condition_to_physical_group(bc, pg_outer_bound)

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

            geo.create_mesh(2)
            mesh = geo.get_mesh(2)
            regions = geo.get_regions()

        # Defining the shape function for solving the FE problem
        prb = Problem("Power Cable", mesh, None, regions, materials, boundary_conditions, excitations)

        # Error: AttributeError: 'TriCartesianEdgeShapeFunction' object has no attribute 'elem2regi'
        '''prb_elec = ElectricProblemCartStatic(self.model_name, shape_function, mesh, regions, materials,
                                             boundary_conditions, excitations)'''


        return prb

# Task 2:  magnetic flux line patterns

def main():

    f = 1e3
    omega = 2 * np.pi * f
    l = 1
    R = np.eye(3) * 1e-3
    G = np.zeros((3, 3))
    L = np.array([[5, 1, 1], [1, 5, 1], [1, 1, 5]]) * 1e-6
    C = np.array([[5, -1, -1], [-1, 5, -1], [-1, -1, 5]]) * 1e-9

    Z = R + 1j * omega * L
    Y = G + 1j * omega * C

    _, Tu = la.eig(Z @ Y)
    _, Ti = la.eig(Y @ Z)

    Zm = la.inv(Tu) @ Z @ Ti
    Ym = la.inv(Ti) @ Y @ Tu

    Ti_real = np.real(Ti)
    Ti_imag = np.imag(Ti)

    print('Tu', Tu[:, 0], Tu[:, 1], Tu[:, 2])
    print('Ti', Ti[:, 0], Ti[:, 1], Ti[:, 2])

    bcs = []
    k = list(np.arange(3))

    for i in range(Ti.shape[0]):

        power_cable = PowerCable(current=Ti_real[:, i])
        problem = power_cable.create_problem(mesh_size_factor=0.2, magn=True, show_gui=False) # bcs
        mesh: TriMesh = problem.mesh
        shape_function = TriCartesianEdgeShapeFunction(mesh)
        problem.shape_function = shape_function

        # ValueError: Matrix A is singular, because it contains empty row(s)
        curlcurl = shape_function.curlcurl_operator(problem.regions, problem.materials, Reluctivity)
        mass = shape_function.mass_matrix(problem.regions, problem.materials, Conductivity)

        matrix = curlcurl + 1j * omega * mass
        load1 = shape_function.load_vector(problem.regions, problem.excitations)
        for e_i, exci in enumerate(problem.excitations):
            exci.value = Ti_imag[e_i, i]

        # da komplexer Strom
        load2 = shape_function.load_vector(problem.regions, problem.excitations)
        load = load1 + 1j * load2

        matrix_shrink, rhs_shrink, _, _, support_data = shape_function.shrink(matrix, load, problem, 1)
        a_shrink, _ = type(problem).solve_linear_system(matrix_shrink.tocsr(), rhs_shrink.tocsr())
        vector_potential = shape_function.inflate(a_shrink, problem, support_data)
        real_vp = np.real(vector_potential)

        b_field = shape_function.curl(real_vp)

        # Plot 1: Wire1 doppelt so großen Wert wie die anderen und nicht komplex
        # Plot 2: alle Wire gleichen Realteil, Wire1 nicht komplex, anderen e-17j
        # Plot 3: Wire sehr sehr gering, Wire 2, 3 fast gleich groß
        mesh.plot_scalar_field(real_vp, title="Vector Potential distr")
        mesh.plot_equilines(real_vp, title='Äquipotetiallinien Vector Potential')
        # nur das hier war die Aufgabe
        mesh.plot_scalar_field(np.linalg.norm(b_field, axis=1), title="Absolute b field")



        plt.show()

if __name__ == '__main__':
    main()
