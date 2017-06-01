# from . basic import *
# from . wf2 import omega0_k, g
from scipy.integrate import quad
import numpy as np

kcutoff = 1e3

#


def nu(gBB):
    return np.sqrt(gBB)


def ur(mI, mB):
    return (mB * mI) / (mB + mI)


def eB(k, mB):
    return k**2 / (2 * mB)


def w(k, gBB, mB, n0):
    return np.sqrt(eB(k, mB) * (eB(k, mB) + 2 * gBB * n0))


def Wk(k, gBB, mB, n0):
    return np.sqrt(eB(k, mB) / w(k, gBB, mB, n0))


def omega0_k(k, gBB, mI, mB, n0):
    return w(k, gBB, mB, n0) + (k**2 / (2 * mI))


def g(aIBi, kcutoff, gBB, mI, mB, n0):
    # assuming finite vector of k values, gives interaction strength constant
    # divSum = (4 * np.pi / (2 * np.pi)**3) * sum((kVec**2) * (Wk(kVec, gBB, mB, n0)**2) / omega_k(P, PB, kVec, gBB, mI, mB, n0))
    return 1 / ((ur(mI, mB) / (2 * np.pi)) * aIBi - (ur(mI, mB) / np.pi**2) * kcutoff)


#

def Bk0(k, aIBi, gBB, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    integrand = lambda kp: (Wk(kp, gBB, mB, n0) / omega0_k(kp, gBB, mI, mB, n0)) * kp**2
    val, abserr = quad(integrand, 0, kcutoff)
    s = (4 * np.pi / (2 * np.pi)**3) * val
    return -1 * np.sqrt(n0) * (Wk(k, gBB, mB, n0) / omega0_k(k, gBB, mI, mB, n0)) * 1 / (g(aIBi, kcutoff, gBB, mI, mB, n0)**(-1) + s)


def num_phonons0(aIBi, gBB, mI, mB, n0):
    integrand = lambda k: (np.abs(Bk0(k, aIBi, gBB, mI, mB, n0))**2) * (k**2)
    val, abserr = quad(integrand, 0, kcutoff)
    return (4 * np.pi / (2 * np.pi)**3) * val


gBB = 0.05
mI = 1
mB = 1
n0 = 1
aIBi = -100

Nph0 = num_phonons0(aIBi, gBB, mI, mB, n0)
print(Nph0)
