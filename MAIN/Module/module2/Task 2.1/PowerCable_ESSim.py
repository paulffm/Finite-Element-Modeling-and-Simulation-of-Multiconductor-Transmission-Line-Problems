import numpy as np
import matplotlib.pyplot as plt
from PowerCable import PowerCable_Elec
from pyrit.bdrycond import BCDirichlet, BdryCond, BCFloating
from pyrit.region import Regions
from utils import get_bc_idx
show_plot = True


'''
Task 2.1 c): Simulation von PowerCable in E-Statik zur Berechnung von C mit low frequency approach
'''

def main():
    ''' HDG 2010AC bei f = 1.0 kHz Low Frequency Approx noch gemacht s.10'''
    phi_elec = []
    K_lst = []
    Xu_lst = []
    u_val = 40
    n = 3

    for i in range(n):
        voltages = np.zeros(n)
        voltages[i] = u_val
        print(voltages)
        bcs = [BCDirichlet(val) if val != 0 else BCFloating() for val in voltages]

        power_cable = PowerCable_Elec()
        problem, shape_function = power_cable.create_problem(bcs, mesh_size_factor=0.2, show_gui=False)

        solution = problem.solve()
        mesh = problem.mesh
        divgrad_matrix = solution.divgrad_matrix
        #print('divgrad', divgrad_matrix.shape)

        #print('Energy', solution.energy)
        K_lst.append(np.asarray(divgrad_matrix.toarray()))
        phi_elec.append(np.asarray(solution.potential).reshape(-1, 1))
        #print('pot', solution.potential.shape)


        # get indices with boundary conditions to shrink phi
        regions_of_bc = problem.boundary_conditions.regions_of_bc(problem.regions)
        '''Dict with the key being the ID of a boundary condition and the value being e list of IDs of the regions that
            have this boundary condition.'''
        '''keyslst = list(regions_of_bc.keys())
        print(keyslst, keyslst[2])

        # floating bc
        idx_floating = shape_function._calc_values_for_floating([keyslst[2]], mesh, problem.regions, regions_of_bc)
        idx_floating2 = shape_function._calc_values_for_floating([keyslst[3]], mesh,  problem.regions, regions_of_bc)
        lst_bc = idx_floating[0].tolist()
        lst_bc.extend(idx_floating2[0])
        print('float', lst_bc)

        # dirchlet idx
        idx_dir, _ = shape_function._calc_values_for_dirichlet(problem.boundary_conditions, [keyslst[0]],
                                                            mesh, problem.regions, regions_of_bc)
        idx_dir2, _ = shape_function._calc_values_for_dirichlet(problem.boundary_conditions, [keyslst[1]],
                                                            mesh, problem.regions, regions_of_bc)
        print('idx123', idx_dir[0], idx_dir2[0])
        lst_dir = idx_dir[0].tolist()
        lst_dir.extend(idx_dir2[0].tolist())
        print('dir', lst_dir)
        lst_bc.extend(lst_dir)
        print('lst_bc', lst_bc)'''

        # get indices with boundary conditions to shrink phi
        lst_bc = get_bc_idx(i)
        phi_dof = np.copy(solution.potential)
        for bc_i in (lst_bc):
            phi_dof[bc_i] = 0
        Xu = (solution.potential - phi_dof) / u_val
        Xu_lst.append(Xu.reshape(-1, 1))

        # Plots the magnetic flux density, style options are 'arrows', 'abs', 'stream'.
        if show_plot:
            solution.plot_d_field('abs')
            #solution.plot_d_field('arrows')
            #solution.plot_d_field('stream')
            solution.plot_equilines()
            solution.plot_potential()
            solution.plot_e_field('abs')
            solution.plot_energy_density()
            plt.show()


    phi_arr = np.hstack((phi_elec[0], phi_elec[1], phi_elec[2]))
    Xu_arr = np.hstack((Xu_lst[0], Xu_lst[1], Xu_lst[2]))

    print('C:', phi_arr.T @ K_lst[0] @ phi_arr / u_val ** 2)
    print('C nach FE-Skript:', Xu_arr.T @ K_lst[0] @ phi_arr / u_val)

    # alle K gleich
    #print(K_lst[0] == K_lst[1])
    '''
    C [[1.08478611e-09 3.13061093e-10 3.13033335e-10]
        [3.13061093e-10 1.08478930e-09 3.12990823e-10]
        [3.13033335e-10 3.12990823e-10 1.08477673e-09]]
    '''

if __name__ == '__main__':
    main()