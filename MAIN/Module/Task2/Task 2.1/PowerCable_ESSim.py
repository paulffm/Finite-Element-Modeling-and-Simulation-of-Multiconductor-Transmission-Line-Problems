import numpy as np
import matplotlib.pyplot as plt
from PowerCable import PowerCable_Elec
from pyrit.bdrycond import BCDirichlet, BdryCond, BCFloating
show_plot = True

'''
Task 2.1 c): Simulation von PowerCable in E-Statik zur Berechnung von C mit low frequency approach
'''

def main():

    phi_elec = []
    K_list = []

    n = 3

    for i in range(n):
        voltages = np.zeros(n)
        voltages[i] = 1
        print(voltages)
        bcs = [BCDirichlet(val) if val != 0 else BCFloating() for val in voltages]

        power_cable = PowerCable_Elec()
        problem = power_cable.create_problem(bcs, mesh_size_factor=0.2, show_gui=False)

        solution = problem.solve()
        mesh = problem.mesh
        divgrad_matrix = solution.divgrad_matrix
        #print('divgrad', divgrad_matrix.shape)

        #print('Energy', solution.energy)
        K_list.append(np.asarray(divgrad_matrix.toarray()))
        phi_elec.append(np.asarray(solution.potential).reshape(-1, 1))
        #print('pot', solution.potential.shape)

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

    Ke_arr = np.hstack((K_list[0], K_list[1], K_list[2]))
    phi_arr = np.hstack((phi_elec[0], phi_elec[1], phi_elec[2]))

    print('C', phi_arr.T @ K_list[0] @ phi_arr)
    '''
    C [[1.08478611e-09 3.13061093e-10 3.13033335e-10]
        [3.13061093e-10 1.08478930e-09 3.12990823e-10]
        [3.13033335e-10 3.12990823e-10 1.08477673e-09]]
    '''

if __name__ == '__main__':
    main()