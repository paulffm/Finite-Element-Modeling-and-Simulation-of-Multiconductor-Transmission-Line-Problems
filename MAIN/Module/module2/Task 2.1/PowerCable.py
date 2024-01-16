# coding=utf-8
from dataclasses import dataclass
import numpy as np
from pyrit.geometry import Geometry, Circle, Surface
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.constants import mu_0, epsilon_0
import gmsh
from pyrit.geometry import Geometry
from pyrit.bdrycond import BCDirichlet, BdryCond, BCFloating
from pyrit.material import Mat, Materials, Reluctivity
from pyrit.material import Materials, Mat, Permittivity, Conductivity, Reluctivity
from pyrit.shapefunction import TriCartesianEdgeShapeFunction, TriCartesianNodalShapeFunction
from pyrit.excitation import Excitations, CurrentDensity
from pyrit.problem import MagneticProblemCartStatic, ElectricProblemCartStatic
from pyrit.toolbox.PostprocessToolbox import plot_field_on_line, get_field_on_line
show_plot = True
"""
Implementation of the class PowerCable that represents a power cable with 3 conductors inside for a Electrostatic and
magnetostatic Problem.
"""


@dataclass
class PowerCable_Magn:

    # given
    wire_radius: float = 1.1e-3
    radius_wires_center_points: float = 1.5e-3
    outer_conductor_inner_radius: float = 3e-3
    outer_conductor_outer_radius: float = 3.2e-3
    conductivity_copper: float = 57.7e6
    conductivity_surrounding: float = 0
    relative_permittivity_surrounding: float = 10
    id_wire_u = 1
    id_wire_v = 2
    id_wire_w = 3
    id_insulation = 10
    id_outer_conductor = 11
    id_outer_bound = 12
    current: float = 1
    model_name: str = "PowerCable_ys"
    depth: float = 1



    @property
    def ids_wires(self):
        return [self.id_wire_u, self.id_wire_v, self.id_wire_w]

    @property
    def current_density(self):
        """Current density in [A/m^2]"""
        return self.current / (np.pi * self.wire_radius ** 2)


    def create_problem(self, k, type, **kwargs):
        geo = Geometry("Power cable", **kwargs)


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


        excitations = Excitations(exci := CurrentDensity(self.current_density))
        if k == 0:
            geo.add_excitation_to_physical_group(exci, pg_wire_u)
        elif k == 1:
            geo.add_excitation_to_physical_group(exci, pg_wire_v)
        elif k == 2:
            geo.add_excitation_to_physical_group(exci, pg_wire_w)

        boundary_conditions_mag = BdryCond(bc := BCDirichlet(0))

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



        shape_function = TriCartesianEdgeShapeFunction(mesh, self.depth)

        # Setting up the FE problem
        prb = MagneticProblemCartStatic(self.model_name, shape_function, mesh, regions, materials,
                                    boundary_conditions_mag, excitations)



        return prb, shape_function

@dataclass
class PowerCable_Elec:

    '''def __init__(self, k):
        self.k = k'''

    # given
    wire_radius: float = 1.1e-3
    radius_wires_center_points: float = 1.5e-3
    outer_conductor_inner_radius: float = 3e-3
    outer_conductor_outer_radius: float = 3.2e-3
    conductivity_copper: float = 57.7e6
    conductivity_surrounding: float = 0
    relative_permittivity_surrounding: float = 10
    id_wire_u = 1
    id_wire_v = 2
    id_wire_w = 3
    id_insulation = 10
    id_outer_conductor = 11
    id_outer_bound = 12
    voltage: float = 1
    current: float = 1
    model_name: str = "PowerCable_ys"
    depth: float = 1

    @property
    def ids_wires(self):
        return [self.id_wire_u, self.id_wire_v, self.id_wire_w]

    @property
    def current_density(self):
        """Current density in [A/m^2]"""
        return self.current / (np.pi * self.wire_radius ** 2)

    def create_problem(self, bcs, **kwargs):
        geo = Geometry("Power cable", **kwargs)

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

        # self.voltage
        boundary_conditions_elec = BdryCond(bc_outer := BCDirichlet(0), *bcs)
        for bc, pg in zip(bcs, pgs_conductor_bound):
            geo.add_boundary_condition_to_physical_group(bc, pg)

        #geo.add_boundary_condition_to_physical_group(bc_outer, pg_outer_bound)
        geo.add_boundary_condition_to_physical_group(bc_outer, pg_outer_conductor)

        #
        # given
        # %% Building model in gmsh and creating the mesh
        with geo:
            # creating the wires?
            circle_u = Circle(self.radius_wires_center_points, 0, 0, self.wire_radius)
            circle_v = Circle(-0.5 * self.radius_wires_center_points,
                              np.sqrt(3) / 2 * self.radius_wires_center_points,
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

        shape_function = TriCartesianNodalShapeFunction(mesh)

        prb = ElectricProblemCartStatic(self.model_name, tri_mesh=mesh, tri_cart_nodal_shape_function=shape_function,
                                        regions=regions,
                                        materials=materials, boundary_conditions=boundary_conditions_elec)

        '''         # Error: AttributeError: 'TriCartesianEdgeShapeFunction' object has no attribute 'elem2regi'
            shape_function = TriCartesianNodalShapeFunction(mesh)
            #shape_function = TriCartesianEdgeShapeFunction(mesh, self.depth)
            prb = ElectricProblemCartStatic(self.model_name, mesh, shape_function, regions, materials,
                                                 boundary_conditions_elec, excitations)'''

        return prb, shape_function
