from polrabi.wf2 import *
import matplotlib
import matplotlib.pyplot as plt


# # Initialization


matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

gBB = 0.05
mI = 10
mB = 1
n0 = 1


# # Dynamics

w_rot = 0
rfreq = 0.1
aIBi_up = -100
# aIBi_down = 100  # actually assume this is Inf when we set g_down = 0 ?
# P = 0
# PB = 0

kcutoff = 10
tMax = 20
dk = 0.05
# dt = 1 / (10 * kcutoff**2)
dt = 1 / (100 * kcutoff**2)

kVec = np.arange(dk, kcutoff, dk)  # number of phonons seems very sensitive to where we start near k=0
tVec = np.arange(0, tMax, dt)

# initial conditions
a0 = np.array([[0.0],
               [1.0]], dtype=complex)
Bk0 = np.zeros(kVec.size, dtype=complex)

# precomputing things that only depend on k and not t
OmegaK = omega0_k(kVec, gBB, mI, mB, n0)
Wkv = Wk(kVec, gBB, mB, n0)
g_up = g(aIBi_up, kcutoff, gBB, mI, mB, n0)
g_down = 0
gVec = np.array([[g_up],
                 [g_down]], dtype=complex)

print(g_up)

# setting initial beta vector and initializing matrices
Bkt = Bk0
ast = a0

amag = np.zeros((2, tVec.size), dtype=complex)
xp = np.zeros(tVec.size, dtype=complex)
xm = np.zeros(tVec.size, dtype=complex)
NB = np.zeros(tVec.size, dtype=float)


for ind, t in enumerate(tVec):

    # keep track of quantities we care about (storing data)

    # print(ast) -- starts getting big for t>16, also overflows for dk<1
    xpt = pchi(Bkt, kVec, gBB, mB, n0)
    xmt = mchi(Bkt, kVec, gBB, mB, n0)
    atemp = np.abs(ast)**2

    amag[[0, 1], ind] = atemp.reshape((2,))
    xp[ind] = xpt
    xm[ind] = xmt
    NB[ind] = np.dot(np.conjugate(Bkt), Bkt)

    # calculate some useful quantities that will be useful later in the loop

    agSum = np.dot(atemp.reshape(2,), gVec.reshape(2,))
    asti = np.array([[ast[1, 0]],
                     [ast[0, 0]]], dtype=complex)
    OmBSum = (4 * np.pi / (2 * np.pi)**3) * np.dot((kVec**2) * OmegaK, np.abs(Bkt)**2)
    wtemp = w_rot * ast[0, 0]
    EsDiff = np.array([[wtemp],
                       [0.0]], dtype=complex)

    # update Bkt and ast to the t+1 value

    BDiff = -1j * (OmegaK * Bkt + (np.sqrt(n0) * Wkv + Wkv * xpt + Wkv**(-1) * xmt) * agSum)
    aDiff = -1j * (ast * gVec * (np.sqrt(n0) + xpt + xmt) * (np.sqrt(n0) + xpt - xmt) - ast * (OmBSum + (np.sqrt(n0) * xpt + xpt**2 - xmt**2) * agSum) + rfreq * asti + EsDiff)

    Bkt = Bkt + dt * BDiff
    ast = ast + dt * aDiff


figS, axS = plt.subplots()
axS.plot(tVec, amag[0, :], 'b-', label=r'$P_{\uparrow}$')
axS.plot(tVec, amag[1, :], 'r-', label=r'$P_{\downarrow}$')
axS.set_xlabel('Time ($t$)')
axS.set_ylabel('Probability')
axS.set_title(r'Probability of being in spin state')
axS.legend()

figC, axC = plt.subplots()
axC.plot(tVec, xp, 'b-', label=r'$\chi^{+}$')
axC.plot(tVec, xm, 'r-', label=r'$\chi^{-}$')
axC.set_xlabel('Time ($t$)')
axC.set_ylabel(r'$\chi(\beta)$')
axC.set_title(r'$\chi(\beta)$')
axC.legend()

figN, axN = plt.subplots()
axN.plot(tVec, NB, 'k-')
axN.set_xlabel('Time ($t$)')
axN.set_ylabel('$N_{ph}$')
axN.set_title('Number of Phonons')

plt.show()
