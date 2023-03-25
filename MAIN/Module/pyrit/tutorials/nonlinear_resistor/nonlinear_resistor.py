# -*- coding: utf-8 -*-
r"""
Problem Description

A cylindrical resistor that is located between two electrodes to which a voltage of 15kV is applied.
The resistor's material has a strongly field-dependent conductivity which has field grading effect.
The geometry is imported from a .geo file with predefined physical groups.

.. sectionauthor:: ruppert
"""
# pylint: disable=unused-import,redefined-outer-name

from dataclasses import dataclass
from typing import Union, TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from pyrit.material import Mat, Materials, Conductivity, DifferentialConductivity
from pyrit.bdrycond import BdryCond, BCDirichlet
from pyrit.toolbox.ImportGmshToolbox import geo_to_msh, read_msh_file, generate_regions
from pyrit.mesh import AxiMesh
from pyrit.shapefunction import TriAxisymmetricNodalShapeFunction
from pyrit.problem import CurrentFlowProblemAxiStatic, CurrentFlowSolverInfoAxiStatic
from pyrit.toolbox.NonlinearMaterialToolbox import NonlinearProperty, NonlinearDifferentialProperty
from pyrit.toolbox.PostprocessToolbox import plot_field_on_line

if TYPE_CHECKING:
    from pathlib import Path


@dataclass
class NonlinearResistor:
    """Example for a stationary current problem with nonlinear materials."""

    # Parameters describing the nonlinear conductivity
    p0: float = 1e-10
    p1: float = 1864
    p2: float = 0.7e6
    p3: float = 2.4e6

    # Excitation parameters
    ground_potential: float = 0  # [V]
    excitation_voltage: float = 15e3  # [V]

    # Parameters for geometry import from gmsh
    path_to_geo_file: Union[str, 'Path'] = 'nonlinear_resistor.geo'

    def sigma(self, solution, element):
        """Nonlinear conductivity function."""
        e_field_abs = solution.e_field_abs[element]
        return self.p0 * ((1 + self.p1 ** ((e_field_abs - self.p2) / self.p2)) /
                              (1 + self.p1 ** ((e_field_abs - self.p3) / self.p2)))

    def differential_sigma(self, solution, element):
        """Derivative dsigma/d(e_field_abs^2)."""
        e_field_abs = solution.e_field_abs[element]
        dfdE = self.p0 * np.log(self.p1) / self.p2 * (self.p1 ** (e_field_abs / self.p2) *
                                                          (1 / self.p1 - self.p1 ** (- self.p3 / self.p2)) /
                                                          (1 + self.p1 ** ((e_field_abs - self.p3) / self.p2)) ** 2)
        return np.divide(dfdE, 2 * e_field_abs, out=np.zeros_like(dfdE), where=e_field_abs != 0)

    def create_problem(self, refinement_steps: int = 0, show_gui: bool = False) -> CurrentFlowProblemAxiStatic:
        """Create the problem.

        Parameters
        ----------
        refinement_steps : int
            Number of refinement steps in the meshing process.
        show_gui: bool
            Open the GMSH GUI?

        Returns
        -------
        problem : CurrentFlowProblemAxiStatic
            A problem instance
        """
        # Materials
        nonlinear_conductivity = NonlinearProperty(self.sigma)
        nonlinear_differential_conductivity = NonlinearDifferentialProperty(self.sigma, self.differential_sigma)
        material = Mat('NONLINEAR_MATERIAL', Conductivity(nonlinear_conductivity),
                       DifferentialConductivity(nonlinear_differential_conductivity))
        materials = Materials(material)

        # Boundary conditions
        bc_ground = BCDirichlet(self.ground_potential, name="GROUND")
        bc_voltage = BCDirichlet(self.excitation_voltage, name="VOLTAGE")
        boundary_cond = BdryCond(bc_ground, bc_voltage)

        # Mesh
        path_to_msh_file = geo_to_msh(self.path_to_geo_file, refinement_steps=refinement_steps, show_gui=show_gui)
        mesh = read_msh_file(path_to_msh_file, mesh_type=AxiMesh)

        # Regions
        regions = generate_regions(path_to_msh_file, materials, boundary_cond)

        # Shape functions
        shape_function = TriAxisymmetricNodalShapeFunction(mesh)

        # Problem
        problem = CurrentFlowProblemAxiStatic("Cylindrical Resistor", mesh, shape_function, regions,
                                              materials, boundary_cond)
        return problem


if __name__ == "__main__":
    # Create problem
    nonlinear_resistor = NonlinearResistor()
    problem = nonlinear_resistor.create_problem(refinement_steps=1, show_gui=True)

    # Solve problem
    # solver_settings = CurrentFlowSolverInfoAxiStatic()
    # solver_settings.tolerance_newton = 1e-12
    # solution = problem.solve(solver_info=solver_settings)
    solution = problem.solve()

    # Plot potential
    solution.plot_potential(kind='abs')

    # Plot field in radial direction
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_field_on_line(problem.mesh, [0.01, 0], [0.03, 0],
                       solution.e_field_abs,
                       lower_bound=0.01, upper_bound=0.03,
                       fig=fig, ax=ax, linestyle='', marker='.')
    plt.show()
