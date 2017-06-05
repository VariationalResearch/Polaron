from . basic import *
from . wf2 import omega0_k, g

kcutoff = 20


def pchi(beta_Vec, Wk_Vec, dV_Vec, gBB, mB, n0):
    # takes finite vector of Beta_{k,theta} values and returns (1/2)*Sum_{k,theta}[W_{k}*(Beta_{k,theta}+Beta_{k,theta}^{*})]
    betaSum = beta_Vec + np.conjugate(beta_Vec)
    return (1 / 2) * np.dot(Wk_Vec, betaSum * dV_Vec)


def mchi(beta_Vec, Wki_Vec, dV_Vec, gBB, mB, n0):
    # takes finite vector of Beta_{k,theta} values and returns (1/2)*Sum_{k,theta}[W_{k}^{-1}*(Beta_{k,theta}-Beta_{k,theta}^{*})]
    betaDiff = beta_Vec - np.conjugate(beta_Vec)
    return (1 / 2) * np.dot(Wki_Vec, betaDiff * dV_Vec)


def PB(beta_Vec, kcos_Vec, dV_Vec, gBB, mB, n0):
    # takes finite vector of Beta_{k,theta} values and returns Sum_{k,theta}[k*cos(theta)*|Beta_{k}|^2]
    betaMag = beta_Vec * np.conjugate(beta_Vec)
    return np.dot(kcos_Vec, betaMag * dV_Vec)


def dynOverlap(NB_vec, phi_Vec):
    # dynamical overlap/Ramsey interferometry signal
    exparg = -1j * phi_Vec - (1 / 2) * NB_vec
    return np.exp(exparg)


def spectFunc(S_Vec, t_Vec):
    # spectral function (Fourier Transform of dynamical overlap)
    tstep = t_Vec[1] - t_Vec[0]
    N = t_Vec.size
    tdecay = 3
    decayFactor = np.exp(-1 * t_Vec / tdecay)
    # decayFactor = 1
    sf = 2 * np.real(np.fft.ifft(S_Vec * decayFactor))
    omega = 2 * np.pi * np.fft.fftfreq(N, d=tstep)
    return omega, sf


# def qPCrit(aIBi, gBB, mI, mB, n0):
#     return mI * nu(gBB) + qPB(mI * nu(gBB), aIBi, gBB, mI, mB, n0)


# def qaSi(DP, gBB, mI, mB, n0):
#     integrand = lambda k: (4 * ur(mI, mB) / (k**2) - ((Wk(k, gBB, mB, n0)**2) / (DP * k / mI)) * np.log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)))) * (k**2)
#     val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
#     return (1 / (2 * np.pi * ur(mI, mB))) * val


# def qPB(DP, aIBi, gBB, mI, mB, n0):
#     integrand = lambda k: ((2 * (w(k, gBB, mB, n0) + (k**2) / (2 * mI)) * (DP * k / mI) + (w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * np.log((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) / (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)))) / ((w(k, gBB, mB, n0) + (k**2) / (2 * mI) - (DP * k / mI)) * (w(k, gBB, mB, n0) + (k**2) / (2 * mI) + (DP * k / mI)) * (DP * k / mI)**2)) * (Wk(k, gBB, mB, n0)**2) * (k**3)
#     val, abserr = quad(integrand, 0, upcutoff, epsabs=0, epsrel=1.49e-12)
#     return n0 / (ur(mI, mB)**2 * (aIBi - qaSi(DP, gBB, mI, mB, n0))**2) * val
