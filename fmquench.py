from polrabi.quench import *
import matplotlib
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from scipy.integrate import trapz


# # Initialization

start = timer()

matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

mI = 1
mB = 1
n0 = 1
gBB = (4 * np.pi / mB) * 0.05

P = 0.5


# # Dynamics


aIBi = 2

kcutoff = 10
dk = 0.05

Ntheta = 50
dtheta = np.pi / (Ntheta - 1)

tMax = 5
dt = 1e-5

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
g = g(aIBi, kcutoff, gBB, mI, mB, n0)
thetaones = np.ones(thetaVec.size)

print(g)

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

    BDiff = -1j * (g * np.sqrt(n0) * Wk_Vec + Bkt * (Omega0K_Vec - kcos_Vec * (P - PB_Vec[ind]) / mI) + g * (Wk_Vec * xpt + Wki_Vec * xmt))
    phiDiff = g * n0 + g * np.sqrt(n0) * xpt + (P**2 - PB_Vec[ind]**2) / (2 * mI)
    Bkt = Bkt + dt * BDiff
    phit = phit + dt * phiDiff

    # print([PBt, xpt, xmt])


S_Vec = dynOverlap(NB_Vec, phi_Vec)
freq_Vec, A_Vec = spectFunc(S_Vec, tVec)

end = timer()

print(end - start)

print(trapz(A_Vec, freq_Vec))

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

fig4, ax4 = plt.subplots()
ax4.plot(tVec, np.abs(S_Vec), 'k-')
ax4.set_xlabel('Time ($t$)')
ax4.set_ylabel(r'$S(t)$')
ax4.set_title('Dynamical Overlap')
fig4.savefig('quench_DynamicalOverlap.pdf')

fig5, ax5 = plt.subplots()
ax5.plot(freq_Vec, A_Vec, 'k-')
ax5.set_xlabel(r'Frequency ($\omega$)')
ax5.set_ylabel(r'$A(\omega)$')
ax5.set_title(r'Spectral Function')
fig5.savefig('quench_SpectralFunction.pdf')

plt.show()
