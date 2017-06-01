from numpy import *
from scipy.integrate import quad, dblquad
import matplotlib.pyplot as plt

# should make a module with general functions and then submodules with
# functions only relevant to specific cases (finite momentum, 1st trial
# wf, etc.)


def nu(gBB):
    return sqrt(gBB / 2)


def ur(mI, mB):
    return (mB * mI) / (mB + mI)


def eB(k, mB):
    return k**2 / (2 * mB)


def w(k, gBB, mB, n0):
    return sqrt(eB(k, mB) * (eB(k, mB) + gBB * n0))


def Wk(k, gBB, mB, n0):
    return sqrt(eB(k, mB) / w(k, gBB, mB, n0))


def PMax(aIBi, gBB, mI, mB, n0):
    return nu(gBB) + PB(nu(gBB), aIBi, gBB, mI, mB, n0)


def Eup(P=0, PB=0, aIBi, aSi, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return ((P**2 + PB**2) / (2 * mI)) + 2 * pi * n0 / (ur(mI, mB) * (aIBi - aSi))


def Edown(P=0, PB=0, aIBi, aSi, w_rot, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return Eup(P, PB, aIBi, aSi, mI, mB, n0) + w_rot


def rMass(P, PB, mI):
    return mI * P / (P - PB)


def Bk(k, aIBi, aSi, gBB, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return -2 * pi * sqrt(n0) * Wk(k, gBB, mB, n0) / (ur(mI, mB) * (aIBi - aSi) * (w(k, gBB, n0) + (k**2) / (2 * mI)))


def aSi0(gBB, mI, mB, n0):
    integrand = lambda k: (2 * ur(mI, mB) / k**2 - (Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI))) * (k**2)
    return (1 / (pi * ur(mI, mB))) * quad(integrand, 0, Inf)


def aSi(DP, gBB, mI, mB, n0):
    integrand = lambda k, theta: (2 * ur(mI, mB) / (k**2) - (Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI) * cos(theta))) * (k**2) * sin(theta)
    return (1 / (2 * pi * ur(mI, mB))) * dblquad(integrand, 0, Inf, lambda: 0, lambda: pi)


def PB(DP, aIBi, gBB, mI, mB, n0):
    integrand = lambda k, theta: ((Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI) * cos(theta))**2) * (k**3) * sin(theta) * cos(theta)
    return n0 / (ur(mI, mB)**2 * (aIBi - aSi(DP, gBB, mI, mB, n0))**2) * dblquad(integrand, 0, Inf, lambda: 0, lambda: pi)


# def aSi1(DP, gBB, mI, mB, n0):
#     integrand = lambda k: (4 * ur(mI, mB) / (k**2) - ((Wk(k, gBB, mB, n0)**2) / (DP * k / mI)) * log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)))) * (k**2)
#     return (1 / (2 * pi * ur(mI, mB))) * quad(integrand, 0, Inf)


# def PB1(DP, aIBi, gBB, mI, mB, n0):
#     integrand = lambda k: ((2 * (w(k, gBB, mB, n0) + (k**2) / (2 * mI)) * (DP * k / mI) + (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)))) / ((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * (DP * k / mI)**2)) * (Wk(k, gBB, mB, n0)**2) * (k**3)
#     return n0 / (ur(mI, mB)**2 * (aIBi - aSi1(DP, gBB, mI, mB, n0))**2) * quad(integrand, 0, Inf)


def DP(DPi, P, aIBi, gBB, mI, mB, n0):
    # define module variable Err,limit
    DPg = DPi
    DPn = DPg + 0.5
    lim = limit

    while abs(DPn - DPg) > Err:
        if lim == 0:
            print('Loop convergence limit reached')
            break

        DPn = abs(P - PB(DPg, aIBi, gBB, mI, mB, n0))
        DPg = DPn
        lim = lim - 1

    return DPn


def cs_overlap(aIBi_up, aIBi_down, gBB, mI, mB, n0):
    aSi0_v = aSi0(gBB, mI, mB, n0)
    integrand = lambda k: (abs(Bk(k, aIBi_up, aSi0_v, gBB, mI, mB, n0))**2 + abs(Bk(k, aIBi_down, aSi0_v, gBB, mI, mB, n0))**2 - 2 * Bk(k, aIBi_up, aSi0_v, gBB, mI, mB, n0) * Bk(k, aIBi_down, aSi0_v, gBB, mI, mB, n0)) * (k**2)
    exparg = -1 * (1 / 2) * (4 * pi / (2 * pi) ** 3) * quad(integrand, 0, Inf)
    return exp(exparg)


def qp_residue(aIBi, aSi, mI, gBB):
    integrand = lambda k: (abs(Bk(k, aIBi, aSi, gBB, mI, mB, n0))**2) * (k**2)
    exparg = -1 * (1 / 2) * (4 * pi / (2 * pi)**3) * quad(integrand, 0, Inf)
    return exp(exparg)


def Hspin(aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0):
    xi = cs_overlap(aIBi_up, aIBi_down, gBB, mI, mB, n0)
    aSi0_v = aSi0(gBB, mI, mB, n0)
    H = array([[Eup(0, 0, aIBi_up, aSi0_v, mI, mB, n0), Omega * xi], [Omega * xi, Edown0(0, 0, aIBi_down, aSi0_v, w_rot, mI, mB, n0)]])
    return matrix(H)
