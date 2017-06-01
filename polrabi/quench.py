from . basic import *
from . wf2 import omega0_k, g


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
