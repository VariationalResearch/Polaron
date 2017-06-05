from polrabi.quench import *
import matplotlib
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import os
from scipy.integrate import trapz
from polrabi.staticfm import PCrit


# # Initialization

matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

# mI = 1
# mB = 1
# n0 = 1
# gBB = (4 * np.pi / mB) * 0.05

# P = 0.5


# # # Dynamics


# aIBi = 2

# kcutoff = 10
# dk = 0.05

# Ntheta = 50
# dtheta = np.pi / (Ntheta - 1)

# tMax = 10
# dt = 1e-5


def dynamics(cParams, gParams, sParams):
    # takes parameters, performs dynamics, and outputs desired observables
    [P, aIBi] = cParams
    [kcutoff, dk, Ntheta, dtheta, tMax, dt] = gParams
    [mI, mB, n0, gBB] = sParams

    kVec = np.arange(dk, kcutoff, dk)
    # thetaVec = np.arange(0, np.pi + dtheta, dtheta)
    thetaVec = np.arange(dtheta, np.pi, dtheta)
    tVec = np.arange(0, tMax, dt)

    # initial conditions
    Bk0_mat = np.zeros((thetaVec.size, kVec.size), dtype=complex)
    Bk0_V = Bk0_mat.reshape(thetaVec.size * kVec.size)
    phi0 = 0 + 0j

    # precomputing things that only depend on k,theta and not t
    Omega0K = omega0_k(kVec, gBB, mI, mB, n0)
    Wkv = Wk(kVec, gBB, mB, n0)
    gnum = g(aIBi, kcutoff, gBB, mI, mB, n0)
    thetaones = np.ones(thetaVec.size)

    Omega0K_mat = np.outer(thetaones, Omega0K)
    Wk_mat = np.outer(thetaones, Wkv)
    dV_mat = (2 * np.pi / (2 * np.pi)**3) * np.outer(dtheta * np.sin(thetaVec), dk * kVec**2)
    kcos_mat = np.outer(np.cos(thetaVec), kVec)

    Omega0K_Vec = Omega0K_mat.reshape(thetaVec.size * kVec.size)
    Wk_Vec = Wk_mat.reshape(thetaVec.size * kVec.size)
    Wki_Vec = 1 / Wk_Vec
    dV_Vec = dV_mat.reshape(thetaVec.size * kVec.size)
    kcos_Vec = kcos_mat.reshape(thetaVec.size * kVec.size)

    # calculate differential equation

    # setting initial beta vector and initializing matrices
    Bkt = Bk0_V
    phit = phi0

    PB_Vec = np.zeros(tVec.size, dtype=float)
    phi_Vec = np.zeros(tVec.size, dtype=complex)
    NB_Vec = np.zeros(tVec.size, dtype=float)

    for ind, t in enumerate(tVec):
        # keep track of quantities we care about (storing data)

        PBt = PB(Bkt, kcos_Vec, dV_Vec, gBB, mB, n0)
        PB_Vec[ind] = PBt
        # print(PBt)
        phi_Vec[ind] = phit

        NBt = np.dot(Bkt * np.conjugate(Bkt), dV_Vec)
        NB_Vec[ind] = NBt
        # print(NBt)

        # calculate some useful quantities that will be useful later in the loop

        xpt = pchi(Bkt, Wk_Vec, dV_Vec, gBB, mB, n0)
        xmt = mchi(Bkt, Wki_Vec, dV_Vec, gBB, mB, n0)

        # update Bkt and ast to the t+1 value

        BDiff = -1j * (gnum * np.sqrt(n0) * Wk_Vec + Bkt * (Omega0K_Vec - kcos_Vec * (P - PB_Vec[ind]) / mI) + gnum * (Wk_Vec * xpt + Wki_Vec * xmt))
        phiDiff = gnum * n0 + gnum * np.sqrt(n0) * xpt + (P**2 - PB_Vec[ind]**2) / (2 * mI)
        Bkt = Bkt + dt * BDiff
        phit = phit + dt * phiDiff

        # print([PBt, xpt, xmt])

    S_Vec = dynOverlap(NB_Vec, phi_Vec)
    freqVec, A_Vec = spectFunc(S_Vec, tVec)

    # save data
    tfData = [tVec, freqVec]
    paramData = [cParams, gParams, sParams]
    obData = [PB_Vec, NB_Vec, S_Vec, A_Vec]
    data = [paramData, tfData, obData]

    dirpath = os.path.dirname(os.path.realpath(__file__))
    np.save(dirpath + '/data/fmquench_aIBi:%.2f_P:%.2f.npy' % (aIBi, P), data)


# calculate dynamics

mI = 1
mB = 1
n0 = 1
gBB = (4 * np.pi / mB) * 0.05

P = 0.85
aIBi = -20
print(PCrit(aIBi, gBB, mI, mB, n0))

kcutoff = 20
dk = 0.05

Ntheta = 10
dtheta = np.pi / (Ntheta - 1)

tMax = 3
dt = 1e-5

cParams = [P, aIBi]
gParams = [kcutoff, dk, Ntheta, dtheta, tMax, dt]
sParams = [mI, mB, n0, gBB]


start = timer()

dynamics(cParams, gParams, sParams)

end = timer()

print(end - start)

# print(trapz(A_Vec, freq_Vec))

# figN, axN = plt.subplots()
# axN.plot(tVec, NB_Vec, 'k-')
# axN.set_xlabel('Time ($t$)')
# axN.set_ylabel('$N_{ph}$')
# axN.set_title('Number of Phonons')
# figN.savefig('quench_PhononNumber.pdf')

# figPB, axPB = plt.subplots()
# axPB.plot(tVec, PB_Vec, 'b-')
# axPB.set_xlabel('Time ($t$)')
# axPB.set_ylabel('$P_{B}$')
# axPB.set_title('Phonon Momentum')
# figPB.savefig('quench_PhononMomentum.pdf')

# figp, axp = plt.subplots()
# axp.plot(tVec, np.sign(phi_Vec) * np.remainder(np.abs(phi_Vec), 2 * np.pi) / np.pi, 'r-')
# axp.set_xlabel('Time ($t$)')
# axp.set_ylabel(r'$\frac{\phi(t)}{\pi}$')
# axp.set_title('Global Phase')
# figp.savefig('quench_GlobalPhase.pdf')

# fig, axes = plt.subplots(nrows=1, ncols=2)
# axes[0].plot(tVec, np.abs(S_Vec), 'k-')
# axes[0].set_xlabel('Time ($t$)')
# axes[0].set_ylabel(r'$\left|S(t)\right|$')
# axes[0].set_title('Dynamical Overlap')


# axes[1].plot(freqVec, A_Vec, 'k-')
# axes[1].set_xlim([-30, 30])
# axes[1].set_ylim([0, 0.1])
# axes[1].set_xlabel(r'Frequency ($\omega$)')
# axes[1].set_ylabel(r'$A(\omega)$')
# axes[1].set_title(r'Spectral Function')
# fig.savefig('quench_DynOverlap&SpectFunction.pdf')

# plt.show()
