from . basic import *
from scipy.integrate import quad
# from scipy.integrate import dblquad
from scipy import interpolate

err = 1e-5
limit = 1e3
alpha = 0.005
upcutoff = 1e3


def PCrit(aIBi, gBB, mI, mB, n0):
    return mI * nu(gBB) + PB(mI * nu(gBB), aIBi, gBB, mI, mB, n0)


def Eup(P, PB, aIBi, aSi, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return ((P**2 - PB**2) / (2 * mI)) + 2 * np.pi * n0 / (ur(mI, mB) * (aIBi - aSi))


def Edown(P, PB, aIBi, aSi, w_rot, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return Eup(P, PB, aIBi, aSi, mI, mB, n0) + w_rot


def rMass(P, PB, mI):
    m = mI * P / (P - PB)

    if np.isscalar(P):
        if P == 0:
            return 1
        else:
            return m
    else:
        mask = (P == 0)
        m[mask] = 1
        return m


def num_phonons(aIBi, aSi, gBB, mI, mB, n0):
    integrand = lambda k: (np.abs(Bk(k, aIBi, aSi, gBB, mI, mB, n0))**2) * (k**2)
    val, abserr = quad(integrand, 0, upcutoff)
    return (4 * np.pi / (2 * np.pi)**3) * val


def qp_residue(aIBi, aSi, gBB, mI, mB, n0):
    exparg = -1 * (1 / 2) * num_phonons(aIBi, aSi, gBB, mI, mB, n0)
    return np.exp(exparg)


def Bk(k, aIBi, aSi, gBB, mI, mB, n0):
    # plug in aSi=aSin0(mI,gBB) for the P=0 case
    return -2 * np.pi * np.sqrt(n0) * Wk(k, gBB, mB, n0) / (ur(mI, mB) * (aIBi - aSi) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI)))


def aSi0(gBB, mI, mB, n0):
    integrand = lambda k: (2 * ur(mI, mB) / k**2 - (Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI))) * (k**2)
    val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
    return (1 / (np.pi * ur(mI, mB))) * val


# def aSi2(DP, gBB, mI, mB, n0):
#     integrand = lambda k, theta: (2 * ur(mI, mB) / (k**2) - (Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI) * cos(theta))) * (k**2) * sin(theta)
#     val, abserr = dblquad(integrand, 0, upcutoff, lambda k: 0, lambda k: np.pi)
#     return (1 / (2 * np.pi * ur(mI, mB))) * val


# def PB2(DP, aIBi, gBB, mI, mB, n0):
#     integrand = lambda k, theta: ((Wk(k, gBB, mB, n0)**2) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI) * cos(theta))**2) * (k**3) * sin(theta) * cos(theta)
#     val, abserr = dblquad(integrand, 0, upcutoff, lambda k: 0, lambda k: np.pi)
#     return n0 / (ur(mI, mB)**2 * (aIBi - aSi2(DP, gBB, mI, mB, n0))**2) * val


def aSi(DP, gBB, mI, mB, n0):
    integrand = lambda k: (4 * ur(mI, mB) / (k**2) - ((Wk(k, gBB, mB, n0)**2) / (DP * k / mI)) * np.log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)))) * (k**2)
    val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
    return (1 / (2 * np.pi * ur(mI, mB))) * val


def PB(DP, aIBi, gBB, mI, mB, n0):
    integrand = lambda k: ((2 * (w(k, gBB, mB, n0) + (k**2) / (2 * mI)) * (DP * k / mI) + (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * np.log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)))) / ((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * (DP * k / mI)**2)) * (Wk(k, gBB, mB, n0)**2) * (k**3)
    val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
    return n0 / (ur(mI, mB)**2 * (aIBi - aSi(DP, gBB, mI, mB, n0))**2) * val


def DP(DPi, P, aIBi, gBB, mI, mB, n0):
    global err, limit, alpha
    DP_old = DPi
    DP_new = 0
    lim = np.copy(limit)

    while True:

        if lim == 0:
            print('Loop convergence limit reached')
            return -1

        DP_new = DP_old * (1 - alpha) + alpha * np.abs(P - PB(DP_old, aIBi, gBB, mI, mB, n0))
        # print(DP_old, DP_new)

        if np.abs(DP_new - DP_old) < err:
            break
        else:
            DP_old = np.copy(DP_new)

        lim = lim - 1

    return DP_new


# interpolation method functions

def PB_integral(DP, gBB, mI, mB, n0):
    integrand = lambda k: ((2 * (w(k, gBB, mB, n0) + (k**2) / (2 * mI)) * (DP * k / mI) + (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * np.log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)))) / ((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * (DP * k / mI)**2)) * (Wk(k, gBB, mB, n0)**2) * (k**3)
    val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
    return val


def PB_prefactor(aIBi, aSi, gBB, mI, mB, n0):
    return n0 / (ur(mI, mB)**2 * (aIBi - aSi)**2)


def createSpline(Nsteps, gBB, mI, mB, n0):
    DP_max = mI * nu(gBB)
    DP_step = DP_max / Nsteps
    DPVals = np.arange(DP_step, DP_max, DP_step)
    aSiVals = np.zeros(DPVals.size)
    PBintVals = np.zeros(DPVals.size)

    for idp, DP in enumerate(DPVals):
        aSiVals[idp] = aSi(DP, gBB, mI, mB, n0)
        PBintVals[idp] = PB_integral(DP, gBB, mI, mB, n0)

    DPVals = np.concatenate((np.array([0]), DPVals))
    aSiVals = np.concatenate((np.array([aSi0(gBB, mI, mB, n0)]), aSiVals))
    PBintVals = np.concatenate((np.array([0]), PBintVals))

    aSi_tck = interpolate.splrep(DPVals, aSiVals, s=0)
    PBint_tck = interpolate.splrep(DPVals, PBintVals, s=0)

    np.save('aSi_spline.npy', aSi_tck)
    np.save('PBint_spline.npy', PBint_tck)


def aSi_interp(DP, aSi_tck):
    return interpolate.splev(DP, aSi_tck, der=0)


def PB_interp(DP, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck):
    aSi = aSi_interp(DP, aSi_tck)
    pf = PB_prefactor(aIBi, aSi, gBB, mI, mB, n0)
    return pf * interpolate.splev(DP, PBint_tck, der=0)


def DP_interp(DPi, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck):
    global err, limit, alpha
    DP_old = DPi
    DP_new = 0
    lim = np.copy(limit)

    while True:

        if lim == 0:
            print('Loop convergence limit reached')
            return -1

        DP_new = DP_old * (1 - alpha) + alpha * np.abs(P - PB_interp(DP_old, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck))
        # print(DP_old, DP_new)

        if np.abs(DP_new - DP_old) < err:
            break
        else:
            DP_old = np.copy(DP_new)

        lim = lim - 1

    return DP_new
