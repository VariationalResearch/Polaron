from . basic import *

# NOTE: FORGOT TO MULTIPLY BY dk IN INTEGRALS FOR SECONDWF CALCULATIONS
# ALSO, WK**-1 IS WRONG APPARENTLY -> NEED TO DO 1/WK


def pchi(betaVec, kVec, gBB, mB, n0):
    # takes finite vector of Beta_{k} values and returns (1/2)*Sum_{k}[W_{k}*(Beta_{k}+Beta_{k}^{*})]
    betaSum = betaVec + np.conjugate(betaVec)
    return (1 / 2) * (4 * np.pi / (2 * np.pi)**3) * np.dot(Wk(kVec, gBB, mB, n0), betaSum * kVec**2)


def mchi(betaVec, kVec, gBB, mB, n0):
    # takes finite vector of Beta_{k} values and returns (1/2)*Sum_{k}[W_{k}^{-1}*(Beta_{k}-Beta_{k}^{*})]
    betaDiff = betaVec - np.conjugate(betaVec)
    return (1 / 2) * (4 * np.pi / (2 * np.pi)**3) * np.dot(Wk(kVec, gBB, mB, n0)**(-1), betaDiff * kVec**2)


# def omega_k(P, PB, k, gBB, mI, mB, n0):
#     return w(k, gBB, mB, n0) + (k**2 / (2 * mI)) - (k / mI) * (P - PB)


def omega0_k(k, gBB, mI, mB, n0):
    return w(k, gBB, mB, n0) + (k**2 / (2 * mI))


def g(aIBi, kcutoff, gBB, mI, mB, n0):
    # assuming finite vector of k values, gives interaction strength constant
    # divSum = (4 * np.pi / (2 * np.pi)**3) * sum((kVec**2) * (Wk(kVec, gBB, mB, n0)**2) / omega_k(P, PB, kVec, gBB, mI, mB, n0))
    return 1 / ((ur(mI, mB) / (2 * np.pi)) * aIBi - (ur(mI, mB) / np.pi**2) * kcutoff)
