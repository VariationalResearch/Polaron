from polrabi.staticfm import *
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
# from timeit import default_timer as timer


# # INITIALIZATION


matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

gBB = 0.05
mI = 1
mB = 1
n0 = 1

nuV = nu(gBB)
DPMax = mI * nuV
res = aSi0(gBB, mI, mB, n0)


# # ZERO MOMENTUM


# aIBi_vals = np.linspace(-1, 1, 1000)
# aSi_v = aSi0(gBB, mI, mB, n0)
# E_vals = Eup(0, 0, aIBi_vals, aSi_v, mI, mB, n0)

# # compare with MATLAB values -- note: biggest discrepancy is near resonance where values are off by < 0.1 (understandable)
# Edat = np.genfromtxt('zm.dat', delimiter=',')
# mask = np.abs(Edat - E_vals) > 1e-2
# print(aIBi_vals[mask])

# fig, ax = plt.subplots()
# ax.plot(aIBi_vals, E_vals, 'k-')
# ax.set_ylim([-50, 50])
# ax.set_xlabel('Inverse Scattering Length ($a_{IB}^{-1}$)')
# ax.set_ylabel('Energy')
# ax.set_title('Polaron Energy at Zero Momentum (static case - no rabi drive)')
# plt.show()


# # INTERPOLATION


# Nsteps = 1e3
# createSpline(Nsteps, gBB, mI, mB, n0)

aSi_tck = np.load('aSi_spline.npy')
PBint_tck = np.load('PBint_spline.npy')

# DP_max = mI * nuV
# DPv = np.linspace(0, DP_max, 100)
# fig, ax = plt.subplots()
# ax.plot(DPv, aSi_interp(DPv, aSi_tck), 'k-')
# ax.plot(DPv, PB_interp(DPv, -4, gBB, mI, mB, n0, aSi_tck, PBint_tck), 'b-')
# plt.show()


# # DP GRID


# aIBiVals = np.linspace(-10, 10, 100)
# PcVals = PCrit(aIBiVals, gBB, mI, mB, n0)

# grid = []

# for ind, aIBi in enumerate(aIBiVals):
#     step = 0.1 * mI * nuV
#     PVals = np.arange(0.1 * mI * nuV, 0.95 * PcVals[ind], step)
#     grid.append((aIBi, PVals))

# points = []
# DPg = 0

# for ind, vertsec in enumerate(grid):
#     (aIBi, PVals) = vertsec
#     for P in PVals:
#         DP_stable = DP_interp(DPg, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#         if (DP_stable == -1 or DP_stable > DPMax):
#             print([aIBi, P, DP_stable])
#             break
#         points.append([aIBi, P, DP_stable])


# pointsarray = np.array(points)
# np.savetxt("aIB_P_DP_points.csv", pointsarray)


# # DATA PROCESSING


# points = np.genfromtxt('aIB_P_DP_points.csv', delimiter=' ')

# aIBip = points[:, 0]
# Pp = points[:, 1]
# DPp = points[:, 2]

# aSip = aSi_interp(DPp, aSi_tck)
# PBp = PB_interp(DPp, aIBip, gBB, mI, mB, n0, aSi_tck, PBint_tck)
# Ep = Eup(Pp, PBp, aIBip, aSip, mI, mB, n0)
# rMp = rMass(Pp, PBp, mI)


# # aIBi, P, DP, aSi, PB, E, rM
# dat = np.concatenate((points, aSip[:, np.newaxis], PBp[:, np.newaxis], Ep[:, np.newaxis], rMp[:, np.newaxis]), axis=1)
# np.savetxt("fmdat.csv", dat)


# # IMPURITY MOMENTUM VS. INTERACTIONS


# aIBiVals = np.linspace(-10, 10, 100)
# Pc_min = PCrit(np.amin(aIBiVals), gBB, mI, mB, n0)
# PVals = np.linspace(0.1 * mI * nuV, 0.95 * Pc_min, 4)

# fig, ax = plt.subplots()
# colortyp = np.array(['r', 'g', 'b', 'y', 'c', 'm', 'y', 'k'])
# ax.plot(aIBiVals, np.zeros(aIBiVals.size), 'k', label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%d$' % 0)

# DPg = 0
# for indp, P in enumerate(PVals):
#     DPVals = np.zeros(aIBiVals.size)
#     for inda, aIBi in enumerate(aIBiVals):
#         DP_stable = DP_interp(DPg, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#         if (DP_stable == -1 or DP_stable > DPMax):
#             print([aIBi, P, DP_stable])
#             DPVals[inda] = float('nan')
#         else:
#             DPVals[inda] = DP_stable
#     mask = ~np.isnan(DPVals)
#     DPValsC = DPVals[mask]
#     aIBiValsC = aIBiVals[mask]

#     DP_tck = interpolate.splrep(aIBiValsC, DPValsC, s=0)
#     DP_int = interpolate.splev(aIBiVals, DP_tck, der=0)
#     Pnorm = P / Pc_min
#     ax.plot(aIBiVals, DP_int / DPMax, colortyp[indp], label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%.2f$' % Pnorm)

# ax.legend()
# ax.set_xlabel(r'Scattering Length ($a_{IB}^{-1}$)')
# ax.set_ylabel(r'Impurity Momentum ($\frac{\Delta P}{m_{I}\nu_{s}}$)')
# ax.set_title('Impurity Momentum vs Interactions')
# plt.show()
# # fig.savefig('impuritymom.pdf')


# # EFFECTIVE MASS VS. INTERACTIONS


aIBiVals = np.linspace(-10, 10, 100)
P = 0.1 * mI * nuV

fig, ax = plt.subplots()
DPg = 0
PBVals = np.zeros(aIBiVals.size)
for inda, aIBi in enumerate(aIBiVals):
    DP_stable = DP_interp(DPg, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
    if (DP_stable == -1 or DP_stable > DPMax):
        print([aIBi, P, DP_stable])
        PBVals[inda] = float('nan')
    else:
        PBVals[inda] = PB_interp(DP_stable, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
mask = ~np.isnan(PBVals)
PBValsC = PBVals[mask]
aIBiValsC = aIBiVals[mask]

ax.plot(aIBiVals, rMass(P, PBVals, mI) / mI, 'b', label=r'$P=0.1 m_{I}\nu_{s}$')
ax.legend()
ax.set_xlabel(r'Scattering Length ($a_{IB}^{-1}$)')
ax.set_ylabel(r'Mass ($\frac{M_{pol}}{m_{I}}=\frac{P}{P-P_{B}}$)')
ax.set_title('Effective Mass vs Interactions')
plt.show()


# # NUMBER OF EXCITATIONS, Z-FACTOR, ENERGY VS INTERACTION STRENGTH

# aIBiVals = np.linspace(-10, 10, 100)
# Pc_min = PCrit(np.amin(aIBiVals), gBB, mI, mB, n0)
# PVals = np.linspace(0.1 * mI * nuV, 0.95 * Pc_min, 4)

# fig, ax = plt.subplots()
# fig2, ax2 = plt.subplots()
# fig3, ax3 = plt.subplots()
# colortyp = np.array(['r', 'g', 'b', 'y', 'c', 'm', 'y', 'k'])
# # ax.plot(aIBiVals, np.zeros(aIBiVals.size), 'k', label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%d$' % 0)

# DPg = 0
# for indp, P in enumerate(PVals):
#     NBVals = np.zeros(aIBiVals.size)
#     qpVals = np.zeros(aIBiVals.size)
#     EVals = np.zeros(aIBiVals.size)
#     for inda, aIBi in enumerate(aIBiVals):
#         DP_stable = DP_interp(DPg, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#         if (DP_stable == -1 or DP_stable > DPMax):
#             print([aIBi, P, DP_stable])
#             qpVals[inda] = float('nan')
#             NBVals[inda] = float('nan')
#             EVals[inda] = float('nan')
#         else:
#             aSi = aSi_interp(DP_stable, aSi_tck)
#             PB = PB_interp(DP_stable, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#             NBVals[inda] = num_phonons(aIBi, aSi, gBB, mI, mB, n0)
#             qpVals[inda] = qp_residue(aIBi, aSi, gBB, mI, mB, n0)
#             EVals[inda] = Eup(P, PB, aIBi, aSi, mI, mB, n0)
#     mask = ~np.isnan(NBVals)
#     NBValsC = NBVals[mask]
#     qpValsC = qpVals[mask]
#     aIBiValsC = aIBiVals[mask]
#     EValsC = EVals[mask]

#     qp_tck = interpolate.splrep(aIBiValsC, qpValsC, s=0)
#     qp_int = interpolate.splev(aIBiVals, qp_tck, der=0)
#     Pnorm = P / Pc_min
#     ax.plot(aIBiVals, NBVals, colortyp[indp], label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%.2f$' % Pnorm)
#     ax2.plot(aIBiVals, qp_int, colortyp[indp], label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%.2f$' % Pnorm)
#     ax3.plot(aIBiVals, EVals, colortyp[indp], label=r'$\frac{P}{P_{crit}(a_{IB}^{-1}=-10)}=%.2f$' % Pnorm)

# ax.legend()
# ax.set_xlabel(r'Scattering Length ($a_{IB}^{-1}$)')
# ax.set_ylabel(r'Number of Phonons ($N_{ph}$)')
# ax.set_title('Number of Phonons vs Interactions')

# ax2.legend()
# ax2.set_xlabel(r'Scattering Length ($a_{IB}^{-1}$)')
# ax2.set_ylabel(r'Quasiparticle Residue ($e^{-\frac{1}{2}N_{ph}}$)')
# ax2.set_title('Quasiparticle Residue vs Interactions')

# ax3.legend()
# ax3.set_xlabel(r'Scattering Length ($a_{IB}^{-1}$)')
# ax3.set_ylabel(r'Energy)')
# ax3.set_title('Energy vs Interactions')

# plt.show()
# fig.savefig('impuritymom.pdf')


# # ENERGY VS MOMENTUM


# aIBiVals = np.array([res + 0.3, res + 0.5, 1, 3, 5, 10])
# # aIBiVals = np.array([-5, res - 0.3, res + 0.3, 5])

# fig, ax = plt.subplots()
# fign, axn = plt.subplots()
# colortyp = np.array(['r', 'g', 'b', 'y', 'c', 'm', 'y', 'k'])

# E1 = []
# E2 = []

# DPg = 0
# for inda, aIBi in enumerate(aIBiVals):
#     Pc = PCrit(aIBi, gBB, mI, mB, n0)
#     PVals = np.linspace(0.1 * mI * nuV, 0.95 * Pc, 100)
#     EVals = np.zeros(PVals.size)
#     for indp, P in enumerate(PVals):
#         DP_stable = DP_interp(DPg, P, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#         if (DP_stable == -1 or DP_stable > DPMax):
#             print([aIBi, P, DP_stable])
#             EVals[indp] = float('nan')
#         else:
#             aSi = aSi_interp(DP_stable, aSi_tck)
#             PB = PB_interp(DP_stable, aIBi, gBB, mI, mB, n0, aSi_tck, PBint_tck)
#             EVals[indp] = Eup(P, PB, aIBi, aSi, mI, mB, n0)
#     mask = ~np.isnan(EVals)
#     EValsC = EVals[mask]
#     PValsC = PVals[mask]

#     a0 = aSi0(gBB, mI, mB, n0)
#     E0 = Eup(0, 0, aIBi, a0, mI, mB, n0)
#     E_tck = interpolate.splrep(PValsC, EValsC - E0, s=0)
#     E_int = interpolate.splev(PVals, E_tck, der=0)
#     E1_int = interpolate.splev(PVals, E_tck, der=1)
#     E2_int = interpolate.splev(PVals, E_tck, der=2)
#     E1.append((PVals / Pc, E1_int))
#     E2.append((PVals / Pc, E2_int))

#     ax.plot(PVals, EVals - E0, colortyp[inda], label=r'$a_{IB}^{-1}=%.2f$' % aIBi)
#     axn.plot(PVals / Pc, EVals - E0, colortyp[inda], label=r'$a_{IB}^{-1}=%.2f$' % aIBi)

# ax.legend()
# ax.set_xlabel('Momentum ($P$)')
# ax.set_ylabel('Energy ($E-E(P=0)$)')
# ax.set_title('Energy vs Momentum')

# axn.legend()
# axn.set_xlabel(r'Momentum ($\frac{P}{P_{crit}(a_{IB})}$)')
# axn.set_ylabel('Energy ($E-E(P=0)$)')
# axn.set_title('Energy vs Momentum')


# fig2, ax2 = plt.subplots()
# (Pnorm, E1Vals) = E1[0]
# (Pnorm, E2Vals) = E2[0]
# ax2.plot(Pnorm, E1Vals, colortyp[0], label=r'$\frac{\partial E}{\partial P}$')
# ax2.plot(Pnorm, E2Vals, colortyp[1], label=r'$\frac{\partial^{2} E}{\partial P^{2}}$')
# ax2.legend()
# ax2.set_xlabel(r'Momentum ($\frac{P}{P_{crit}(a_{IB})}$)')
# ax2.set_ylabel('Energy Derivatives')
# ax2.set_title(r'Energy Behavior for $a_{IB}^{-1}=%.2f$' % aIBiVals[0])


# fig3, ax3 = plt.subplots()
# (Pnorm, E1Vals) = E1[-2]
# (Pnorm, E2Vals) = E2[-2]
# ax3.plot(Pnorm, E1Vals, colortyp[0], label=r'$\frac{\partial E}{\partial P}$')
# ax3.plot(Pnorm, E2Vals, colortyp[1], label=r'$\frac{\partial^{2} E}{\partial P^{2}}$')
# ax3.legend()
# ax3.set_xlabel(r'Momentum ($\frac{P}{P_{crit}(a_{IB})}$)')
# ax3.set_ylabel('Energy Derivatives')
# ax3.set_title(r'Energy Behavior for $a_{IB}^{-1}=%.2f$' % aIBiVals[-2])

# plt.show()
# # fig.savefig('impuritymom.pdf')
