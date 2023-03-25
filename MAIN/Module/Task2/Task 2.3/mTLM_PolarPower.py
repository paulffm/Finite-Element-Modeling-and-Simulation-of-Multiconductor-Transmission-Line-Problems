import numpy as np
from matplotlib import pyplot as plt
from mTLM import A_m, A, solve, P, Pout
from matplotlib.pyplot import cm
import numpy.linalg as la
import scipy.linalg as sla
show_plot = True
'''
Task 2.2a) Konstanten
Task 2.2b) Phasor Diagramme von Spannungen und Ströme
Task 2.3a) Berechnung u0, ul, i0, il und Berechnung von Pin und Pout
Task 2.3b) Berechnung von Power out über Frequenz

'''
def Ak_m(Zkm_ch: float, bkm: float, l: float) -> np.ndarray:
    """Calculates the local modal propagation matrix for the given tlm parameters."""
    bl = bkm * l
    return np.array([
        [np.cosh(bl), -Zkm_ch * np.sinh(bl)],
        [-np.sinh(bl) / Zkm_ch, np.cosh(bl)]
    ])

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

    # (ZY) Tu = lambda Tu
    # (YZ) Ti = lambda Ti
    _, Tu = la.eig(Z @ Y)
    _, Ti = la.eig(Y @ Z)

    # Zm = Tu^-1 * Z * Ti
    # Ym = Ti^-1 * Y * Tu
    Zm = la.inv(Tu) @ Z @ Ti
    Ym = la.inv(Ti) @ Y @ Tu

    # Widerstand: wie breitet sich Mode aus: daher diag
    Zm_ch = np.sqrt(np.diag(Zm) / np.diag(Ym))

    # Ausbreitungskonstante für jeden Mode
    bm = np.diag(np.sqrt(Zm * Ym))

    # Modale Propagationsmatrix:
    ak_m = [Ak_m(Zkm, bkm, l) for Zkm, bkm in zip(Zm_ch, bm)]

    # Am ist Blockdiagonale von ak_m
    am = sla.block_diag(ak_m[0], ak_m[1], ak_m[2])

    # umsortieren sodass u1, u2, u3, i1, i2, i3 in dieser Reihenfolge steht
    idx = [0, 2, 4, 1, 3, 5]
    am = am[idx, :][:, idx]

    # A = Qui * Am Qui^-1
    Q = sla.block_diag(Tu, Ti)
    print(Q.shape)
    print(am.shape)
    a = Q @ am @ la.inv(Q)

    '''print('Z', Z_z, Z_z.shape)   #3x3
    print('Y', Y_y, Y_y.shape)      #3x3
    print('Zm', Zm, Zm.shape)       #3x3
    print('Ym', Ym, Ym.shape)       #3x3
    print('Tu', Tu, Tu.shape)       #3x3
    print('Ti', Ti, Ti.shape)       #3x3
    print('z_char', z_char, z_char.shape) #3x1
    print('beta', b, b.shape)             #3x1
    print('AK_m', Ak_m(z_char[0], b[0], l).shape) # 2x2
    print('am', am, am.shape, am[0], am[0].shape)             #6x6
    print('a', a, a.shape)                #6x6'''


    ## Task 2:

    # a) Ausbreitungskonstante = attenuation const + j phase const = alpha + j beta = np.sqrt(R+jwL) * (G+jwC))
    # alpha = Dämpfungskonstante beta = Phasendrehung

    print('Modes along the cable:')
    print(' Zm:  ', np.diag(Zm))
    print(' Ym:  ', np.diag(Ym))
    print(' beta:', bm)
    print(' Z_char:', Zm_ch)
    print(' phase const:   ', np.imag(bm))
    print(' attenuation const:', np.real(bm))
    print('')


    '''Task 2.2b): Phasor Diagramme Spannungen und Ströme'''
    if show_plot:
        plt.figure()
        color = iter(cm.rainbow(np.linspace(0, 1, 10)))
        for i in range(Tu.shape[1]):
            # um im: (1 0 0) (0 1 0) (0 0 1)
            u_i = Tu[:, i]
            i_i = Ti[:, i]
            for j in range(u_i.shape[0]):
                c = next(color)
                plt.polar([0, np.angle(u_i[j], deg=True)], [0, np.abs(u_i[j])], marker='o', label=f'I{j+1}, U{j+1} Mode {i+1}', c=c)
                plt.polar([0, np.angle(i_i[j], deg=True)], [0, np.abs(i_i[j])], marker='x', c=c)
            plt.legend(loc='lower left').set_draggable(True)
            plt.title(f'Polar Plot of U and I for the 3 Modes: o=Voltage, x=Current')
            plt.show()

            '''for k in range(u_i.shape[0]):
                c = next(color)
                plt.polar([0, np.angle(u_i[k], deg=True)], [0, np.abs(u_i[k])], marker='o',
                          label=f'I{k + 1}, U{k + 1} Mode {i + 1}', c=c)
                plt.polar([0, np.angle(i_i[k], deg=True)], [0, np.abs(i_i[k])], marker='x', c=c)
                plt.legend(loc='lower left').set_draggable(True)
                plt.title(f'Polar Plot of U and I for the 3 Modes: o=Voltage, x=Current')
                plt.show()'''



    ### Task 3 ###

    '''Task 2.3a) Berechnung u0, ul, i0, il und Berechnung von Pin und Pout'''
    l_3 = 2e3
    # Anregung U(0) als Phasor: /np.sqrt(2) und Winkel: cos(wt-phi) => e^-phi
    u0 = np.array([100, 80 * np.exp(np.pi * (-2j / 3)), 60 * np.exp(np.pi * (- 4j / 3))])

    # build System:
    # Bedingungen auf Blatt:
    'M11, M12:'
    # (U(l), I(l)) = A * (U(0), I(0)) => A * (U(0), I(0)) - I_6 * (U(l), I(l)) = 0_6
    'M21:'
    # U(0) =U(0) => I_3 U(0) = U(0)
    'M22:'
    # U(l) = I(l) => I_3 U(l) - I_3 I(l) = 0_6

    zero_33 = np.zeros([3, 3])
    zero_3 = np.array([0, 0, 0])
    # Widerstand R = 1 am Abschluss: Es folgt: U(l) = I(l)
    r = 1

    M11 = a
    M12 = -np.eye(6)
    M21 = sla.block_diag(np.eye(3), zero_33)
    M22 = np.block([[zero_33, zero_33], [np.eye(3), -np.diag([r, r, r])]])
    M = np.block([[M11, M12], [M21, M22]])
    b = np.concatenate([zero_3, zero_3, u0, zero_3])

    # Solve the System
    x = la.solve(M, b)
    u0 = x[0:3]
    i0 = x[3:6]
    ul = x[6:9]
    il = x[9:12]

    p_in = np.real(0.5 * np.sum(u0 * np.conj(i0)))
    p_out = np.real(0.5 * np.sum(ul * np.conj(il)))
    print(f"Power loss {p_in - p_out}W and relative loss {100 * (p_in - p_out) / p_in} %")

    '''Task 2.3b) Berechnung von Power out über Frequenz'''
    # same as before: but in functions
    f_list = np.logspace(1, 6, 200)
    p_out = [Pout(fi, l_3, R, G, L, C, u0) for fi in f_list]

    if show_plot:
        plt.loglog(f_list, p_out)
        plt.xlabel('frequency')
        plt.ylabel('Power out')
        plt.show()


if __name__ == '__main__':
    main()
