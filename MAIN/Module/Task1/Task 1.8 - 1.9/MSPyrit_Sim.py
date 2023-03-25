import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.constants import mu_0
import gmsh
from pyrit.geometry import Geometry
from pyrit.bdrycond import BCDirichlet, BdryCond
from pyrit.material import Mat, Materials, Reluctivity
from pyrit.shapefunction import TriCartesianEdgeShapeFunction
from pyrit.excitation import Excitations, CurrentDensity
from pyrit.problem import MagneticProblemCartStatic
from pyrit.toolbox.PostprocessToolbox import plot_field_on_line, get_field_on_line

# Variable to show or suppress the plots.
show_plot = True


@dataclass
class WireProblem:
    """Class representing the wire problem.

    The geometry is a wire with radius r_1, enclosed within a cylindrical, non-conducting insulating shell with inner
    radius r_1 and outer radius r_2. The outer surface of the insulation shell is considered as a perfect electric
    conductor. The depth of the wire can be given.
    """

    r_1: float = 2e-3  #: inner radius (wire) [m]
    r_2: float = 3.5e-3  #: outer radius (shell) [m]
    depth: float = 300e-3  #: Depth of wire (lz) [m]
    model_name: str = "Wire_ys"  #: model name

    mu_shell: float = 5 * mu_0  #: permeability of shell [H/m]
    mu_wire: float = mu_0  #: permeability of wire [H/m]
    current: float = 16.  #: applied current [A]

    region_wire: int = 2  #: Region ID of the wire
    region_shell: int = 3  #: Region ID of the shell

    @property
    def current_density(self):
        """Current density in [A/m^2]"""
        return self.current / (np.pi * self.r_1 ** 2)

    @property
    def analytic_energy(self):
        """Analytic energy in [J]"""
        return self.current ** 2 * self.depth / (4 * np.pi) * (mu_0 / 4 + self.mu_shell * np.log(self.r_2 / self.r_1))

    @property
    def analytic_inductance(self):
        """Analytic inductance in [Vs/A]"""
        return 2 * self.analytic_energy / self.current ** 2

    def h_analytic(self, radius):
        """Analytic solution of magnetic field in [A/m]

        Parameters
        ----------
        radius : float
            The radius.

        Returns
        -------
        out : float
            The magnetic field strength at the given point.
        """

        def h_inside():
            return self.current_density / 2 * radius

        def h_outside():
            return self.current_density / 2 * self.r_1 ** 2 / radius

        condition = radius < self.r_1
        return condition * h_inside() + (~condition) * h_outside()

    def a_analytic(self, radius):
        """Analytic solution of magnetic vector potential in [Vs/m]

        Parameters
        ----------
        radius : float
            The radius.

        Returns
        -------
        out : float
            The vector potential at the given point.
        """

        def a_inside():
            return -self.current / (2 * np.pi) * (
                    mu_0 / 2 * (radius ** 2 - self.r_1 ** 2) / self.r_1 ** 2 + self.mu_shell * np.log(
                self.r_1 / self.r_2))

        def a_outside():
            return -self.mu_shell * self.current / (2 * np.pi) * np.log(radius / self.r_2)

        condition = radius < self.r_1
        return condition * a_inside() + (~condition) * a_outside()

    def create_problem(self, mesh_size: float = 0.2, refinement_steps: int = 0, show_gui: bool = False):
        """Create the Problem instance

        Parameters
        ----------
        mesh_size : float, optional
            The global mesh size factor for gmsh. Default is 0.2
        refinement_steps : int, optional
            The number of mesh refinement steps. Default is 0
        show_gui : bool, optional
            Boolean that indicates if the gmsh gui should be opened after geometry creation or not. Default is False

        Returns
        -------
        problem : Problems.MagneticProblem2D
            A problem instance
        """
        # %% Initialize instance of gmsh handler
        geo = Geometry(self.model_name, show_gui=show_gui, mesh_size_factor=mesh_size)
        # more options are available to control the mesh
        # using 'mesh_size_factor' to control grid size here.

        # %% Creating materials
        materials = Materials(mat_shell := Mat('shell', Reluctivity(1 / self.mu_shell)),
                              mat_wire := Mat('wire', Reluctivity(1 / self.mu_wire)))

        # %% Creating boundary conditions
        boundary_conditions = BdryCond(bc := BCDirichlet(0))  # Dirichlet-type boundary condition
        # for the magnetic vector potential

        # %% Creating the excitation
        excitations = Excitations(exci := CurrentDensity(self.current_density))

        # %% Creating physical groups
        # A physical group is a set of geometrical entities of the same dimension (see gmsh).
        outer_bound = geo.create_physical_group(1, 1, 'outer_bound')  # grounded outer boundary
        wire = geo.create_physical_group(self.region_wire, 2, 'wire')  # conductive wire carrying the excitation current
        shell = geo.create_physical_group(self.region_shell, 2, 'shell')  # insulating material

        # %% Assigning materials, boundary conditions and excitations to physical groups
        geo.add_material_to_physical_group(mat_shell, shell)
        geo.add_material_to_physical_group(mat_wire, wire)
        geo.add_boundary_condition_to_physical_group(bc, outer_bound)
        geo.add_excitation_to_physical_group(exci, wire)

        # %% Building model in gmsh and creating the mesh.
        with geo:
            # creating the wire with specified radii.
            inner_disc = gmsh.model.occ.addDisk(0, 0, 0, self.r_1, self.r_1)
            outer_disc = gmsh.model.occ.addDisk(0, 0, 0, self.r_2, self.r_2)

            # create the shell; the second entity
            # (wire inner) intersects the shell.
            gmsh.model.occ.cut([[2, outer_disc]], [[2, inner_disc]], removeTool=False)

            # assigning geometry entities to physical groups.
            outer_bound.add_entity(2)
            wire.add_entity(1)
            shell.add_entity(2)

            # # Alternative commands to control mesh parameters and format.
            # gmsh.option.setNumber('Mesh.MeshSizeFactor', 0.1)  # control grid size here.
            # gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)  # MATLAB compatible mesh file format
            # gmsh.model.occ.synchronize()

            # Meshing
            geo.create_mesh(dim=2)
            for _ in range(refinement_steps):
                geo.refine_mesh()
            mesh = geo.get_mesh(dim=2)
            regions = geo.get_regions()

            # if you like to save the created model, use gmsh.write('putname.step')

        # Defining the shape function for solving the FE problem
        shape_function = TriCartesianEdgeShapeFunction(mesh, self.depth)

        # Setting up the FE problem
        prb = MagneticProblemCartStatic(self.model_name, shape_function, mesh, regions, materials,
                                        boundary_conditions, excitations)

        return prb


"""
The problem we want so solve is defined in the class 'WireProblem'. 
With its method 'create_problem' we get an instance of 'MagneticProblem2D', that defines all information needed to
solve the FE problem. This can be done by calling the method 'solve_static' (since we have a static problem). This
method returns a solution object, that again contains all information needed to work with the solution (e.g. to 
calculate energies and fields or to plot something) 
"""

# region Modelling and Simulation
# %% Initialize instance of problem handler
wire_problem = WireProblem()

# %% Modelling and assembly of FE problem
problem = wire_problem.create_problem()

# %% Solving the FE problem
solution = problem.solve()
# endregion

# region Compare energies

# Access variables for post-processing
mesh = problem.mesh
curlcurl_matrix = solution.curlcurl_matrix

# Analytic vector potential on the nodes of the mesh
a_analytic = wire_problem.a_analytic(np.linalg.norm(mesh.node, axis=1))  #: [Tm]

# Magnetic energy with analytic vector potential and numerical curlcurl-matrix
energy_test = wire_problem.depth ** 2 * 1 / 2 * a_analytic @ curlcurl_matrix @ a_analytic  #: [J]
print('L', 2 * solution.energy / wire_problem.current ** 2, 2 * wire_problem.analytic_energy / wire_problem.current ** 2)
# Output important validation QoIs
print('Verification Check:')
print(f'Magnetic energy (analytic solution)                 : {wire_problem.analytic_energy}J')
print(f'Magnetic energy (analytic Az, numerical Knu)        : {energy_test}J')
print(f'Magnetic energy (FE solution)                       : {solution.energy}J')
print('\n')
print('Energy in regions:')
print(f"Energy in region 'wire':  {solution.energy_in_regions(wire_problem.region_wire)}.")
print(f"Energy in region 'shell': {solution.energy_in_regions(wire_problem.region_shell)}.")
print(f"Sum of both:              {solution.energy_in_regions(wire_problem.region_wire, wire_problem.region_shell)}.")

# check the structure of the curlcurl matrix
if show_plot:  # Print structure of matrix
    plt.spy(curlcurl_matrix, markersize=1)
    plt.show()

# endregion

# region Create plots with the provided post-processing routines.

if show_plot:
    # Plots the equilines of the vector potential
    solution.plot_equilines()

    # Plots the magnetic flux density, style options are 'arrows', 'abs', 'stream'.
    solution.plot_b_field('abs')

    # Compute the magnetic flux density magnitude
    b_abs = np.linalg.norm(solution.b_field, axis=1)

    # Plot a field a long an arbitrary line, with/without smoothing option
    start_point = [0, 0]  # x,y
    end_point = [wire_problem.r_2, 0.]  # x,y

    fig, ax = plot_field_on_line(mesh, start_point, end_point, b_abs, lower_bound=0, upper_bound=end_point[0],
                                 smoothing=True, regions=problem.regions, marker='o', label='With smoothing')
    fig, ax = plot_field_on_line(mesh, start_point, end_point, b_abs, lower_bound=0, upper_bound=end_point[0],
                                 smoothing=False, regions=problem.regions, marker='o', fig=fig, ax=ax,
                                 label='Without smoothing')
    ax.legend()
    ax.set_title('Absolute magnetic flux density over the radius')
    ax.set_xlabel('Radius in m')
    ax.set_ylabel(r'$\left|\vec{B}\right|$')
    plt.show()

    # Compare analytic and numerical solution
    radius = np.linspace(1e-5, wire_problem.r_2, 200)
    _, ax = plot_field_on_line(mesh, start_point, end_point, solution.vector_potential_field_values, lower_bound=0,
                               upper_bound=end_point[0], smoothing=True, regions=problem.regions, label='Numeric')
    ax.plot(radius, wire_problem.a_analytic(radius), label='Analytic')
    ax.legend()
    ax.set_title("Vector potential")
    plt.show()

    # if plot_field_on_line() does not fulfill your requirements on the plot you have to create the plot yourself
    # for example if you want to call get_field_on_line() with other default values
    # here is an example on how to use the matplotlib library for plotting

    h_abs = np.linalg.norm(solution.h_field, axis=1)
    _, ax = plt.subplots()

    # extract and interpolate the field along the line.
    points_on_line, normalized_vec_pot = get_field_on_line(mesh, start_point, end_point, 200,
                                                           solution.vector_potential_field_values / max(
                                                               solution.vector_potential_field_values))
    # plot
    x_coordinates = points_on_line[:, 0]
    ax.plot(x_coordinates, normalized_vec_pot, label=r'$A_z$')

    # extract, interpolate incl. smoothing
    points_on_line, normalized_h_field = get_field_on_line(mesh, start_point, end_point, 200, h_abs / max(h_abs),
                                                           smoothing=True, regions=solution.regions)
    # plot
    x_coordinates = points_on_line[:, 0]
    ax.plot(x_coordinates, normalized_h_field, label=r'$H_\varphi$', color=(0.1, 0.5, 0.2))

    ax.set_xlabel("Radius in m")
    ax.set_ylabel(r"Normalized $H_\varphi$ and $A_z$")
    plt.title('Normalized vector potential and magnetic field strength (both numeric)')
    ax.legend()
    plt.show()

    # Pyrit provides selected plotting methods for important fields:
    # 1) post-processing check for magnetic flux density
    solution.plot_b_field('arrows')
    solution.plot_b_field('stream')
    # 2) post-processing check for magnetic field strength
    solution.plot_h_field('abs')
    # 3) post-processing check for magnetic energy
    solution.plot_energy_density()
    # 4) post-processing check for current density
    solution.plot_current_density()
    plt.show()

# region Convergence

"""
We automate the FE procedure such that the calculations are repeated for models with different mesh sizes, 
created by different values for the mesh size option.
"""

refinements_steps = [0, 1, 2, 3]  # variable for the mesh size refinement
energies = np.empty_like(refinements_steps, dtype=float)  # outputs for the convergence plot
element_sizes = np.empty_like(refinements_steps, dtype=float)

"""Iterate over refinement steps. In every step:
- Solve the Problem
- Compute the energy
- Determine the maximum element size
"""
print('Refinement')
for k, refinement in enumerate(refinements_steps):
    problem = wire_problem.create_problem(mesh_size=1, refinement_steps=refinement)
    solution = problem.solve()
    energies[k] = solution.energy
    element_sizes[k] = np.max(solution.mesh.edge_length)
    solution.plot_vector_potential()

# Relative error in the energy and convergence order
rel_error = np.abs((energies - wire_problem.analytic_energy) / wire_problem.analytic_energy)
print(f'rel_errors: {rel_error * 100} %')
print(f'max element sizes: {element_sizes}')
conv_order = (np.log(rel_error[0]) - np.log(rel_error[-1])) / (np.log(element_sizes[0]) - np.log(element_sizes[-1]))
print(f'Convergence order: {conv_order:.2f}')

plt.figure()
plt.loglog(element_sizes, rel_error, label="Error")
plt.xlabel('element size')
plt.ylabel('relative error')
plt.title(f"Relative error. Convergence of order {conv_order:.2f}")
plt.show()

