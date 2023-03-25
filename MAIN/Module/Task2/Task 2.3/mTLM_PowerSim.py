import numpy as np
from matplotlib import pyplot as plt
from mTLM import Z, Y, decompose, Z_ch, beta_m, A_m, A, solve, P, Pout
from matplotlib.pyplot import cm
show_plot = True
'''
Task 2.2a) Konstanten
Task 2.2b) Phasor Diagramme von Spannungen und Ströme
Task 2.3a) Berechnung u0, ul, i0, il und Berechnung von Pin und Pout
Task 2.3b) Berechnung von Power out über Frequenz

'''

def main():

    f = 1e3
    l = 1
    R = np.eye(3) * 1e-3
    G = np.zeros((3, 3))
    L = np.array([[5, 1, 1], [1, 5, 1], [1, 1, 5]]) * 1e-6
    C = np.array([[5, -1, -1], [-1, 5, -1], [-1, -1, 5]]) * 1e-9

    Z_z = Z(R, L, f)
    Y_y = Y(G, C, f)

    Zm, Ym, Tu, Ti = decompose(Z_z, Y_y)

    z_char = Z_ch(Zm, Ym)

    b = beta_m(Zm, Ym)

    am = A_m(z_char, b, l)

    a = A(am, Tu, Ti)


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
    print(' beta:', b)
    print(' Z_char:', z_char)
    print(' phase const:   ', np.imag(b))
    print(' attenuation const:', np.real(b))
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
    l_3 = 2e3  #andere Resonanzen
    u0 = np.array([100, 80 * np.exp(np.pi * (-2j / 3)), 60 * np.exp(np.pi * (- 4j / 3))])
    r = 1
    u0, i0, ul, il = solve(a, u0, r)
    p_in, p_out = P(u0, i0, ul, il)
    print(f"Power loss {p_in - p_out}W and relative loss {100 * (p_in - p_out) / p_in} %")

    '''Task 2.3b) Berechnung von Power out über Frequenz'''
    f_list = np.logspace(1, 6, 200)
    p_out = [Pout(fi, l_3, R, G, L, C, u0) for fi in f_list]
    if show_plot:
        plt.loglog(f_list, p_out)
        plt.xlabel('frequency')
        plt.ylabel('Power out')
        plt.show()


if __name__ == '__main__':
    main()
