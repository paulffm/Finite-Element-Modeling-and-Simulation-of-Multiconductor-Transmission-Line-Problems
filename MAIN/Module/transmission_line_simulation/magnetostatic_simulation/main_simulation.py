import gmsh
from matplotlib import pyplot as plt
from Mesh import Mesh
from shape_function import ShapeFunction_N
import analytic_sol
from scipy.sparse.linalg import spsolve
import numpy as np
from scipy.sparse import csr_matrix
import plot_properties
from plot_properties import entity_in_physical_group
show_plot = True

def Knu_elem(k, shape_function, reluctivity_elems):
    '''
    :param k:
    :param shape_function:
    :param reluctivity_elems:
    :return:
    '''

    # integral(v rot(Wi) * rot(Wj) dV) =>
    # v (Wi dx,
    # sum over elements: ((b.T * b + c.T * c ) / (4 * area * l_z)) * reluctivity
    # produces 3x172 * 172x3 = 3x3 output
    _, b, c, S = shape_function.get_coeff()

    Knu_e = (np.array(b[k, :], ndmin=2).T @ np.array(b[k, :], ndmin=2) +
           np.array(c[k, :], ndmin=2).T @ np.array(c[k, :], ndmin=2)) \
                      / (4*S[k]*shape_function.depth)
    return reluctivity_elems[k] * Knu_e

def calc_bfield(a, shape_function, msh):

    # bx = sum(c * A / 2 * area) / l_z , by = sum(b * A / 2 * area)  / (l_z)
    _, b, c, S = shape_function.get_coeff()
    b_field = np.vstack([np.sum(c * a[msh.elem_to_node[:]] / (2 * S.reshape(-1, 1)), 1)
                         / shape_function.depth,
                         - np.sum(b * a[msh.elem_to_node[:]] / (2 * S.reshape(-1, 1)), 1)
                         / shape_function.depth]).T
    return b_field


def main():

    model_name = "wire"
    r_1 = 2e-3  # [m]
    r_2 = 3.5e-3  # [m]
    depth = 300e-3  # [m]
    mu_0 = 4 * np.pi * 1e-7  # [H/m]
    mu_shell = 5 * mu_0
    eps_0 = 8.8541878128 * 1e-12
    I = 16  # [A]
    J_0 = I / (np.pi * r_1 ** 2)  # [A/m]

    ##### multi_transmission_line_simulation: Construction of a finite-element model with Gmsh #####
    msh = Mesh.create()
    #print(msh.num_elements, msh.num_node): 172, 101
    pg = gmsh.model.getPhysicalGroups()
    # (dim, tag)
    # pg [(1, 3), (2, 1), (2, 2)]

    # Physical Groups
    # In every entry in physical_groups we define the following structure (TAG (dimension, name, indices of all nodes))
    # getNodesforPhysicalGroup(dim, tag)
    physical_groups = dict()

    ##### machine_slot_simulation: Visualization of mesh and regions #####

    for group in pg:
        # physical_groups unterscheidung durch tags: zugriff: physical_groups[3]
        # getphysicalname(dim, tag) -> name as string 'GND'
        # getNodesforPhysicalGroup(dim, tag) -> node_tags, node_coord hier nur node_tags wegen [0];

        physical_groups[group[1]] = (group[0], gmsh.model.getPhysicalName(group[0], group[1]),
                                     gmsh.model.mesh.getNodesForPhysicalGroup(group[0], group[1])[0] - 1)

        # Dict: (Tag, (Dim, physicalName, array with node tags, array with coord of nodes)
        # 3 GND 1 WIRE 2 SHELL
        # access by tag physical_groups[1] oder [2] [3] dann ein ganzes dict
        # physical_groups[1][2] nur tags von wire
        # physical_groups[2][2] name von wire

    #print(physical_groups[1][2])
    gmsh.finalize()

    # indices for all elements by physical group
    elem_shell = entity_in_physical_group(physical_groups, msh.elem_to_node, 'SHELL')
    elem_wire = entity_in_physical_group(physical_groups, msh.elem_to_node, 'WIRE')

    # reluctivity in elem
    # Permeabilität = Durchlässigkeit Magnetfeld; Reluktanz=magn. Widerstand
    reluctivity_elem = 1 / mu_0 * np.ones(msh.num_elements)  # [m/H]
    reluctivity_elem[elem_shell] = 1 / mu_shell  # [m/H] :

    # Task 4: setup the FE shape functions and assemble the stiffness matrix and load vector.
    # construct shape_function

    # Instantiation of Shape function, calculating of a, b, c, S: num_elem x 3 , num_elem x 1
    shape_function = ShapeFunction_N(depth)
    shape_function.calc_coeff(msh)

    # Assign Knu for global indices: exactly taken from supporting remarks
    idx_row = np.zeros((9 * msh.num_elements), dtype='int')
    idx_col = np.zeros((9 * msh.num_elements), dtype='int')
    elem_entries = np.zeros((9 * msh.num_elements))

    # loop over elements
    for k in range(msh.num_elements):

        # idx für nodes von element k: (3, ): Wie zeile
        global_indices = msh.elem_to_node[k, :]

        # (3, 3) -> zeile dreimal untereinander
        triple_global_indices = np.array([global_indices, global_indices, global_indices])

        # np.reshape(triple_global_indices, (9)) schreibt in eine Zeile: (9, )
        # schreibt nach und nach alle 9 indices nebeneinander und zwar wie folgt:
        # immer 3 mal nebeinander: index_row [44 58 49 44 58 49 44 58 49 ...],
        idx_row[9 * k:9 * k + 9] = np.reshape(triple_global_indices, (9))

        # dann col: [44 44 44 58 58 58 49 49 49  ..]
        idx_col[9 * k:9 * k + 9] = np.reshape(triple_global_indices.T, (9))

        # lokales Knu: ((b.T * b + c.T * c ) / (4 * area * l_z)) * reluctivity -> produziert 3x3 output

        # aus 3x3 wird 1x9 und das in elem_entries geschrieben
        elem_entries[9 * k:9 * k + 9] = np.reshape(Knu_elem(k, shape_function, reluctivity_elem), (9))

    # Zuweisung
    idx_row = idx_row.T
    idx_col = idx_col.tolist()
    elem_entries = elem_entries.tolist()

    # Knu Matrix: # [1/H]
    Knu = csr_matrix((elem_entries, (idx_row, idx_col)))
    print('Knu shape', Knu, Knu.shape)


    ##### Task 4: Define load vector j_grid, element-wise current density #####

    # current density in elementen: num_elem x 1
    j_elems = np.zeros(msh.num_elements)  # [A/m^2]:

    # np.sum(shape_function.element_area[elem_in_wire]) surface area vom wire
    # current_density_elems: jedes Wire dreieck gleichen Anteil an I
    j_elems[elem_wire] = I / np.sum(shape_function.element_area[elem_wire])

    'grid currents_e = J integral(Ne,i dA) = J Ae/3 -> Ni aufintegriert über Fläche für jedes Element'
    # num_elems x 1
    grid_currents = j_elems * shape_function.element_area / 3
    x_elems = grid_currents / I

    # vector with values for current contribution of each element on the nodes.
    values = np.zeros(msh.num_node)
    x_values = np.zeros(msh.num_node)

    # Iteration durch jedes Element j
    # addiere Grid current von Element j zu Node i, wenn Node i Teil von Element i ist
    # Am Ende: Für jeden Node stehen dort die addierten Strombeiträge von jedem Element,
    # in welchem sich Node j befindet

    for i in range(3):
        for j in range((msh.num_elements)):
            idx_node = msh.elem_to_node[j, i]
            values[idx_node] += grid_currents[j]
            x_values[idx_node] += x_elems[j]

    # Zuweisung der Werte zu jeweiligen nodes in sparse format: num_nodes x 1
    j_grid = csr_matrix((values, (np.arange(msh.num_node), np.zeros(msh.num_node))), shape=(msh.num_node, 1))
    x_grid = csr_matrix((x_values, (np.arange(msh.num_node), np.zeros(msh.num_node))), shape=(msh.num_node, 1))
    #print('x', x_grid, x_grid * I == j_grid)

    print('unit of Knu: [1/H] : circuit-reluctance matrix')
    print('unit of load vector: [A/m^2]: current density')

    ##### Task 5: First validation is the check of the magnetic energy #####

    x = np.array(msh.node[:, 0]).T
    y = np.array(msh.node[:, 1]).T

    # Radial Koordinate [m]:
    r = np.sqrt(x ** 2 + y ** 2)

    # Analytische Lösung von A in zur Projektion [Tm]:
    A_analytic = depth * analytic_sol.A_z(r)

    # Magnetische Energie (Analytisch und numerisch) in [J]
    # Gesamtenergie eines magnetostatischen Feldes: W = 0.5 Integral (A * J) dV= 0.5 * Integral (A * K * A) dV (KA=J)
    W_magn_test = 1 / 2 * A_analytic @ Knu @ A_analytic
    W_magn_analytic = analytic_sol.W_magn()

    print('Magnetische Energie (analytisch)           :', W_magn_analytic, 'J')
    print('Magnetische Energie (analytisch Az, numerisch Knu)  :', W_magn_test, 'J')

    ##### Task 6: setup and solve the magnetostatic problem #####
    a = np.zeros((msh.num_node, 1))  # Initialize vector of dofs

    # indices of GND are the boundary:
    idx_bc = physical_groups[3][2] # take only the indices out of dict: 28

    value_bc = np.zeros((len(idx_bc), 1))

    # the indices where nothing is given: indices where we have to calculate a: 73 -> num_nodes - dof = index_constraint
    # restlichen indizes -> DoF: An diesen muss a berechnet werden
    idx_dof = np.setdiff1d(np.arange(msh.num_node), idx_bc).tolist()

    # Reduce the system: Knu (num_nodes x num_nodes) --> (num_dof x num_dof): Remove the BC indizes -> known entries
    Knu_red = Knu[idx_dof, :]
    Knu_red = Knu_red[:, idx_dof]

    # reduce rhs
    j_grid_dof = j_grid[idx_dof]
    rhs = j_grid_dof

    # Solve the system: Ka = j spsolver: Ax=b
    a_shrink = spsolve(Knu_red, rhs)

    # Inflate A back to full size: num_nodes x 1
    a_shrink = np.array(a_shrink, ndmin=2).T  # (73, ) only non boundary nodes
    a[idx_dof] = a_shrink  # filled to (101,1)
    a[idx_bc] = value_bc
    a = a.reshape(len(a))  # (101) for every node

    # L: Calculation
    print('Analytical L', analytic_sol.Inductance())
    print('FE L', x_grid.T @ (a / I), a.T @ Knu.toarray() @ a / I ** 2)

    ##### Task 7: Calculate magnetic flux density #####

    # 172 x 2 num_elements
    b_field = calc_bfield(a, shape_function, msh)
    # num_elements x 1
    b_field_abs = np.linalg.norm(b_field, axis=1)  # [T]: magnitude of the magnetic flux density

    # Compare Results: magnetic energy in [J]
    W_magn_fe = 1 / 2 * a @ Knu @ a
    W_magn_fe2 = np.sum(1 / 2 * reluctivity_elem * b_field_abs ** 2
                        * shape_function.element_area * shape_function.depth)
    print('Validierung:')
    print('Magnetic energy (Analytisch)                :', W_magn_analytic, 'J')
    print('Magnetic energy (FE)                      :', W_magn_fe, 'J')
    print('Magnetic energy (Integrierte FE-LSG)   :', W_magn_fe2, 'J')

    # rel error
    rel_error = np.abs((W_magn_fe - W_magn_analytic) / W_magn_analytic)
    print(f'Relative error of energy: {rel_error}')

    if show_plot:

        # plot of analytical solution
        r_list = np.linspace(0, r_2, 100)
        r_list[0] = r_list[1] * 0.5
        fig, axs = plt.subplots(3, 1)
        axs[0].plot(r_list, analytic_sol.B_phi(r_list))
        axs[0].set_title("B")
        axs[0].set(xlabel='r', ylabel='B')
        axs[1].plot(r_list, analytic_sol.H_phi(r_list))
        axs[1].set_title("H")
        axs[1].set(xlabel='r', ylabel='H')
        axs[2].plot(r_list, analytic_sol.A_z(r_list))
        axs[2].set_title("A")
        axs[2].set(xlabel='r', ylabel='A')

        # plot mesh
        plot_properties.plot_mesh(msh)

        # plot regions of mesh
        plot_properties.plot_regions_of_mesh(msh, physical_groups)

        # plot reluctivity
        plot_properties.plot_reluctivity(msh, reluctivity_elem)

        # plot current density in elements
        plot_properties.plot_current(msh, j_elems)

        # Struktur Knu Matrix
        plt.spy(Knu, markersize=1)
        plt.show()

        # plot sol: on ground = 0
        plot_properties.plot_sol(msh, a)

        # plot b_field:
        plot_properties.plot_bfield(msh, b_field_abs)


if __name__ == '__main__':
    main()

