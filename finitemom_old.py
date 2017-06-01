from polrabi.staticfm import *
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
from timeit import default_timer as timer


# # INITIALIZATION


matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

gBB = 0.05
mI = 1
mB = 1
n0 = 1

aIBi = 3

nuV = nu(gBB)


# # ZERO MOMENTUM


# aIBi_vals = linspace(-1, 1, 1000)
# aSi_v = aSi0(gBB, mI, mB, n0)
# E_vals = Eup(0, 0, aIBi_vals, aSi_v, mI, mB, n0)

# # compare with MATLAB values -- note: biggest discrepancy is near resonance where values are off by < 0.1 (understandable)
# Edat = genfromtxt('zm.dat', delimiter=',')
# mask = abs(Edat - E_vals) > 1e-2
# print(aIBi_vals[mask])

# fig, ax = plt.subplots()
# ax.plot(aIBi_vals, E_vals, 'k-')
# ax.set_ylim([-50, 50])
# ax.set_xlabel('Inverse Scattering Length ($a_{IB}^{-1}$)')
# ax.set_ylabel('Energy')
# ax.set_title('Polaron Energy at Zero Momentum (static case - no rabi drive)')
# plt.show()


# # FINITE MOMENTUM


# has trouble converging near P=0 (DP -> 0 so PB integral blows up) and near resonance (similar reason?)
# P = 1e-3  # pb0 fine when P>0.05 (at least far from resonance)
# aIBi = 10
# DPi = 0.1  # can't make this too small (<0.1), sensitive to this near P=0


# a = aSi0(gBB, mI, mB, n0)
# # dp0 = DP(DPi, P, aIBi, gBB, mI, mB, n0)  # depends on PB and hence aSi functions
# # aS0 = aSi(dp0, gBB, mI, mB, n0)  # stabley consistant with analytic formula for P=0
# # pb0 = PB(dp0, aIBi, gBB, mI, mB, n0)  # consistent for aIBi far from resonance, but problems closer?

# aS0, pb0 = aSi_PB(DPi, P, aIBi, gBB, mI, mB, n0)

# # print(dp0)
# print(a, aS0, pb0)


# # MASS AND ENERGY CALCULATIONS

# aIBi = -3
# Pc = PCrit(aIBi, gBB, mI, mB, n0)
# N = 10
# # PVals = linspace(1 / N * Pc, (N - 1) / N * Pc, N)
# PVals = linspace(0.1 * mI * nuV, 0.95 * Pc, N)
# aSiVals = zeros(N)
# PBVals = zeros(N)
# MVals = zeros(N)
# EVals = zeros(N)

# DPi = 0.9 * PVals[0]

# for x in range(N):
#     Pt = PVals[x]

#     DP_stable = DP(DPi, Pt, aIBi, gBB, mI, mB, n0)
#     if DP_stable == -1:
#         break
#     aSit, PBt = aSi(DP_stable, gBB, mI, mB, n0), PB(DP_stable, aIBi, gBB, mI, mB, n0)

#     # aSit, PBt = aSi_PB(DPi, Pt, aIBi, gBB, mI, mB, n0)
#     aSiVals[x] = aSit
#     PBVals[x] = PBt
#     MVals[x] = rMass(Pt, PBt, mI)
#     EVals[x] = Eup(Pt, PBt, aIBi, aSit, mI, mB, n0)

#     DPi = copy(DP_stable)


# a0 = aSi0(gBB, mI, mB, n0)
# E0 = Eup(0, 0, aIBi, a0, mI, mB, n0)

# tck = interpolate.splrep(PVals, EVals, s=0)
# Eder2 = interpolate.splev(PVals, tck, der=2)
# Eder1 = interpolate.splev(PVals, tck, der=1)

# # print(PVals / Pc)
# # print(MVals / mI)
# # print(PVals**2 / (2 * (EVals - E0)))
# # print(Eder)
# # print(Eder1)

# fig, ax = plt.subplots()
# ax.plot(PVals / Pc, MVals / mI, 'k-', label=r'$\frac{m_{I}\cdot P}{P-P_{B}}$')
# ax.plot(PVals / Pc, PVals**2 / (mI * 2 * (EVals - E0)), 'b-', label=r'$\frac{P^{2}}{2(E-E_{0})}$')
# ax.plot(PVals / Pc, 1 / (Eder2 * mI), 'g-', label=r'$(\frac{\partial^{2} E}{\partial P^{2}})^{-1}$')
# # ax.set_ylim([amin(MVals), amax(MVals)])
# ax.legend()
# ax.set_xlabel(r'Polaron Momentum ($\frac{P}{P_{crit}}$)')
# ax.set_ylabel(r'Polaron Mass ($\frac{M}{m_{I}}$)')
# ax.set_title('Different Calculations of Polaron Mass')

# fig2, ax2 = plt.subplots()
# ax2.plot(PVals / Pc, PVals**2 / (2 * MVals[0]), 'k-', label=r'$\frac{P^{2}}{2M_{P \approx 0}$')
# ax2.plot(PVals / Pc, (EVals - E0), 'b-', label=r'$E-E_{0}$')
# ax2.legend()
# ax2.set_xlabel(r'Polaron Momentum ($\frac{P}{P_{crit}}$)')
# ax2.set_ylabel(r'Polaron Energy')
# ax2.set_title('Different Calculations of Polaron Energy')
# plt.show()


# # RESONANCE TRACKING


# aIB_Mlim = -1.6
# aIB_Plim = 2.2

# aIBiValsM = linspace(-4, aIB_Mlim, 10)
# aIBiValsP = linspace(aIB_Plim, 4, 10)
# aIBiVals = concatenate((aIBiValsM, aIBiValsP))

# Pc_min = PCrit(amin(aIBiValsM), gBB, mI, mB, n0)
# print(Pc_min)
# PVals = linspace(0.1 * mI * nuV, 0.95 * Pc_min, 5)

# aS0 = aSi0(gBB, mI, mB, n0)
# aSiVals = zeros((PVals.size, aIBiVals.size))


# DPi = 0.9 * PVals[0]

# for row, P in enumerate(PVals):
#     for col, aIBi in enumerate(aIBiVals):
#         DP_stable = DP(DPi, P, aIBi, gBB, mI, mB, n0)
#         if DP_stable == -1:
#             break
#         aSiVals[row, col] = aSi(DP_stable, gBB, mI, mB, n0)
#         DPi = copy(DP_stable)


# fig3, ax3 = plt.subplots()
# ax3.plot(aIBiVals, aS0 * ones(aIBiVals.size), 'k-', label='$P=0$')

# colortyp = array(['r', 'g', 'b', 'y', 'c', 'm', 'y', 'k'])

# for row, P in enumerate(PVals):
#     ax3.plot(aIBiVals, aSiVals[row, :], colortyp[row], label='$P=%.4f$' % P)

# ax3.legend()
# ax3.set_xlabel(r'Scattering length ($a_{IB}^{-1}$)')
# ax3.set_ylabel(r'Scattering resonance parameter ($a_{*}^{-1}$)')
# ax3.set_title('Tracking Feshback Resonance')
# plt.show()
# fig3.savefig('resonanceshift.pdf')


# # DP GRID (LEFT & RIGHT)


# start = timer()

# aIB_Mlim = 0.27
# aIB_Plim = 0.38


# aIBiValsM = concatenate((linspace(-10, -4, 5), linspace(-4.1, aIB_Mlim, 20)))
# aIBiValsP = concatenate((linspace(10, 4, 5), linspace(4.1, aIB_Plim, 20)))

# PcValsM = PCrit(aIBiValsM, gBB, mI, mB, n0)
# PcValsP = PCrit(aIBiValsP, gBB, mI, mB, n0)

# gridM = []
# gridP = []


# for ind, aIBi in enumerate(aIBiValsM):
#     step = 0.1 * mI * nuV
#     PVals = arange(0.1 * mI * nuV, 0.95 * PcValsM[ind], step)
#     gridM.append((aIBi, PVals))

# for ind, aIBi in enumerate(aIBiValsP):
#     step = 0.1 * mI * nuV
#     PVals = arange(0.1 * mI * nuV, 0.95 * PcValsP[ind], step)
#     gridP.append((aIBi, PVals))

# pointsM = []
# DPg = 0.9 * (0.1 * mI * nuV)
# DPf = 0

# for ind, vertsec in enumerate(gridM):
#     (aIBi, PVals) = vertsec
#     for indP, P in enumerate(PVals):
#         DP_stable = DP(DPg, P, aIBi, gBB, mI, mB, n0)
#         if (DP_stable == -1):
#             print(aIBi, P)
#             break
#         pointsM.append([aIBi, P, DP_stable])
#         DPg = copy(DP_stable)
#         if(indP == 0):
#             DPf = copy(DP_stable)
#         # print(DP_stable)
#     DPg = DPf


# pointsP = []
# DPg = 0.9 * (0.1 * mI * nuV)
# DPf = 0

# for ind, vertsec in enumerate(gridP):
#     (aIBi, PVals) = vertsec
#     for indP, P in enumerate(PVals):
#         DP_stable = DP(DPg, P, aIBi, gBB, mI, mB, n0)
#         if (DP_stable == -1):
#             print(aIBi, P)
#             break
#         pointsP.append([aIBi, P, DP_stable])
#         DPg = copy(DP_stable)
#         if(indP == 0):
#             DPf = copy(DP_stable)
#         # print(DP_stable)
#     DPg = DPf

# pointsP.reverse()

# points = pointsM + pointsP

# pointsarray = array(points)
# savetxt("aIB_P_DP_points.csv", pointsarray)

# end = timer()
# print(end - start)


# # DATA PROCESSING


# points = genfromtxt('aIB_P_DP_points.csv', delimiter=' ')
