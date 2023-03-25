from typing import Tuple
import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
from numpy import pi, sqrt

def Z(R: np.ndarray, L: np.ndarray, f: float) -> np.ndarray:
    """The impedance Z."""
    return R + 2j * pi * f * L


def Y(G: np.ndarray, C: np.ndarray, f: float) -> np.ndarray:
    """The admittance Y."""
    return G + 2j * pi * f * C


def decompose(Z: np.ndarray, Y: np.ndarray) -> [np.ndarray, np.ndarray]:
    """Calculates the modal decomposition of Z and Y."""
    # Tu Eigvector of Z @ Y und Ti of Y @ Z
    u, Tu = la.eig(Z @ Y)
    i, Ti = la.eig(Y @ Z)


    Zm = la.inv(Tu) @ Z @ Ti
    Ym = la.inv(Ti) @ Y @ Tu
    return Zm, Ym, Tu, Ti


def Z_ch(Zm: np.ndarray, Ym: np.ndarray) -> np.ndarray:
    """The characteristic impedance values for each cable for the given modal Z and Y."""
    return sqrt(np.diag(Zm) / np.diag(Ym))


def beta_m(Zm: np.ndarray, Ym: np.ndarray) -> np.ndarray:
    """Calculates the propagation constants for each cable gamma for the given modal Z and Y."""
    return np.diag(sqrt(Zm * Ym))


def Ak_m(Zkm_ch: float, bkm: float, l: float) -> np.ndarray:
    """Calculates the local modal propagation matrix for the given tlm parameters."""
    bl = bkm * l
    return np.array([
        [np.cosh(bl), -Zkm_ch * np.sinh(bl)],
        [-np.sinh(bl) / Zkm_ch, np.cosh(bl)]
    ])


def A_m(Zm_ch: np.ndarray, bm: np.ndarray, l: float) -> np.ndarray:
    """Calculates the modal propagation matrix.
    :param Zm_ch: Array of characteristic impedance values.
    :param bm: Array of propagation constants.
    :param l: Length of the cable.
    """
    ak_m = [Ak_m(Zkm, bkm, l) for Zkm, bkm in zip(Zm_ch, bm)]

    a_m = sla.block_diag(ak_m[0], ak_m[1], ak_m[2])

    # umsortieren sodass u1, u2, u3, i1, i2, i3 in dieser Reihenfolge steht
    idx = [0, 2, 4, 1, 3, 5]

    return a_m[idx, :][:, idx]


def A(A_m: np.ndarray, Tu: np.ndarray, Ti: np.ndarray) -> np.ndarray:
    """Calculates the propagation matrix."""
    T = sla.block_diag(Tu, Ti)
    return T @ A_m @ la.inv(T)

def solve(A: np.ndarray, u0: np.ndarray, R: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solves the mtlm problem.
    :param A: The transmission matrix.
    :param u0: The source voltage vector.
    :param R: The load resistance.
    :return: The solution vectors u0, i0, ul, il as a tuple.
    """
    '''Erklärung: u0 und il gegeben: 
    '''
    Z = np.zeros([3, 3])
    z = np.array([0, 0, 0])
    M11 = A
    M12 = -np.eye(6)
    M21 = sla.block_diag(np.eye(3), Z)
    M22 = np.block([[Z, Z], [np.eye(3), -np.diag([R, R, R])]])
    M = np.block([[M11, M12], [M21, M22]])
    b = np.concatenate([z, z, u0, z])

    x = la.solve(M, b)
    return x[0:3], x[3:6], x[6:9], x[9:12]

def P(u0, i0, ul, il) -> Tuple[float, float]:
    """Calculates the incoming and outgoing power of the cable."""
    # 0.5
    Sin = 0.5 * np.sum(u0 * np.conj(i0))
    Sout = 0.5 * np.sum(ul * np.conj(il))
    return Sin.real, Sout.real


def Pout(f: float, l: float, R: np.ndarray, G: np.ndarray, L: np.ndarray, C: np.ndarray, u0: np.ndarray) -> float:
    """Calculates the load power for the given tlm parameters and the given frequency."""
    z = Z(R, L, f)
    y = Y(G, C, f)
    Zm, Ym, Tu, Ti = decompose(z, y)
    z_char = Z_ch(Zm, Ym)
    b = beta_m(Zm, Ym)
    am = A_m(z_char, b, l)
    a = A(am, Tu, Ti)
    r = 1
    u0, i0, ul, il = solve(a, u0, r)
    _, p_out = P(u0, i0, ul, il)
    return p_out