from . basic import *
from . staticfm import Bk, aSi0, Eup, Edown
from scipy.integrate import quad
from scipy.linalg import expm

# note: integrals are calculated to k=inf, not upper cutoff


def cs_overlap(aIBi_up, aIBi_down, gBB, mI, mB, n0):
    aSi0_v = aSi0(gBB, mI, mB, n0)
    integrand = lambda k: (np.abs(Bk(k, aIBi_up, aSi0_v, gBB, mI, mB, n0))**2 + np.abs(Bk(k, aIBi_down, aSi0_v, gBB, mI, mB, n0))**2 - 2 * Bk(k, aIBi_up, aSi0_v, gBB, mI, mB, n0) * Bk(k, aIBi_down, aSi0_v, gBB, mI, mB, n0)) * (k**2)
    val, abserr = quad(integrand, 0, np.Inf)
    exparg = -1 * (1 / 2) * (4 * np.pi / (2 * np.pi) ** 3) * val
    return np.exp(exparg)



def Hspin(aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0):
    xi = cs_overlap(aIBi_up, aIBi_down, gBB, mI, mB, n0)
    aSi0_v = aSi0(gBB, mI, mB, n0)
    H = np.array([[Eup(0, 0, aIBi_up, aSi0_v, mI, mB, n0), Omega * xi], [Omega * xi, Edown(0, 0, aIBi_down, aSi0_v, w_rot, mI, mB, n0)]])
    return H


def wfcoeff(t, a0, aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0):
    # returns 2x1 vector with the amplitudes of the wavefunction coefficients (a_{up}(t) and a_{down}(t)) at time t given initial condition (vector) a0
    H = Hspin(aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0)
    return np.dot(expm(-1j * H * t), a0)
