import numpy as np

r_1 = 2e-3                # [m]
r_2 = 3.5e-3              # [m]
depth = 300e-3            # [m]
model_name = "wire"

I = 16                    # [A]
J_0 = I/(np.pi*r_1**2)    # [A/m]
mu_0 = 4*np.pi*1e-7       # [H/m]
mu_shell = 5*mu_0


def H_phi(r):
    '''
    Analytic solution of Magnetic Field
    :param r:
    :return:
    '''


    def H_phi_i(r):
        return J_0 / 2 * r

    def H_phi_a(r):
        return J_0 / 2 * r_1 ** 2 / r

    condition = r < r_1
    return condition * H_phi_i(r) + (~condition) * H_phi_a(r)

def B_phi(r):
    '''
    :param r:
    :return:
    '''
    def B_phi_i(r):
        return J_0 * mu_0 / 2 * r

    def B_phi_a(r):
        return J_0 * mu_shell / 2 * r_1 ** 2 / r

    condition = r < r_1
    return condition * B_phi_i(r) + (~condition) * B_phi_a(r)


def A_z(r):
    '''
    Analytic solution of magnetic vector potential
    :param r:
    :return:
    '''

    def A_z_i(r):
        return -I / (2 * np.pi) * (mu_0 / 2 * (r ** 2 - r_1 ** 2) / r_1 ** 2 + mu_shell * np.log(r_1 / r_2))

    def A_z_a(r):
        return -mu_shell * I / (2 * np.pi) * np.log(r / r_2)

    condition = r < r_1
    return condition * A_z_i(r) + (~condition) * A_z_a(r)

def W_magn():
    # magnetic energy(analytical)
    return I ** 2 * depth / (4 * np.pi) * (mu_0 / 4 + mu_shell * np.log(r_2 / r_1))


def Inductance():
    #[H]   : inductance (analytical)
    return 2 * W_magn() / I ** 2
