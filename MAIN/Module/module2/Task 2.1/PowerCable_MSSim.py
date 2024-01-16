import numpy as np
import matplotlib.pyplot as plt
from PowerCable import PowerCable_Magn
show_plot = True
'''
Task 2.1 a, b): Simulation von PowerCable in M-Statik zur Berechnung von R und L mit low frequency approach
'''
def main():
    '''HDG 2010AC bei f = 1.0 kHz Low Frequency Approx noch gemacht'''
    k = list(np.arange(3))
    X_mag = []
    K_list = []
    a_mag = []
    current_list = []

    for i in k:
        power_cable = PowerCable_Magn()

        problem, shape_function = power_cable.create_problem(k=i, type='magn', mesh_size_factor=0.2, show_gui=False)
        load = shape_function.load_vector(problem.regions, problem.excitations)
        X = load / power_cable.current

        # ValueError: Matrix A is singular, because it contains empty row(s) -> without Reluktanz
        solution = problem.solve()
        mesh = problem.mesh
        curlcurl_matrix = solution.curlcurl_matrix

        print('Energy', solution.energy)
        X_mag.append(np.asarray(X.toarray()))
        K_list.append(np.asarray(curlcurl_matrix.toarray()))
        a_mag.append(np.asarray(solution.vector_potential).reshape(-1, 1))
        current_list.append(power_cable.current)

        # Compute the magnetic flux density magnitude
        b_abs = np.linalg.norm(solution.b_field, axis=1)

        # Plots the magnetic flux density, style options are 'arrows', 'abs', 'stream'.
        if show_plot:
            solution.plot_b_field('abs')
            #solution.plot_b_field('arrows')
            #solution.plot_b_field('stream')
            solution.plot_vector_potential()
            solution.plot_equilines()
            solution.plot_current_density()
            solution.plot_energy_density()
            plt.show()

    Xm_arr = np.concatenate((X_mag[0], X_mag[1], X_mag[2]), axis=1)
    am_arr = np.concatenate((a_mag[0], a_mag[1], a_mag[2]), axis=1)


    r_w = 1.1e-3
    sigma = 57.7e6

    R = np.eye(3) * (1 / (sigma * np.pi * r_w ** 2))
    print('Resistance:', R)
    L = Xm_arr.T @ am_arr / current_list[0]
    print('Inductivity', L, am_arr.T @ K_list[0] @ am_arr)

    '''Inductivity [[2.11889983e-07 6.48576562e-08 6.48568935e-08]
                    [6.48576562e-08 2.11890109e-07 6.48566612e-08]
                    [6.48568935e-08 6.48566612e-08 2.11884862e-07]]'''


if __name__ == '__main__':
    main()