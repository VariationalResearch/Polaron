import numpy as np
# if I import polrabi.basic in a script and also import numpy *, are double copies of numpy functions created?


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


# def Wki(k, gBB, mB, n0):
#     return np.sqrt(w(k, gBB, mB, n0) / eB(k, mB))
