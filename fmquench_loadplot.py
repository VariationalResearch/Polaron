import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from polrabi.basic import nu


# # Initialization

matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

# load data

dirpath = os.path.dirname(os.path.realpath(__file__))

P = 0.85
aIBi = -5
data = np.load(dirpath + '/data/fmquench_aIBi:%.2f_P:%.2f.npy' % (aIBi, P))

[paramData, tfData, obData] = data
[cParams, gParams, sParams] = paramData
[tVec, freqVec] = tfData
[PB_Vec, NB_Vec, S_Vec, A_Vec] = obData


# plot

figN, axN = plt.subplots()
axN.plot(tVec, NB_Vec, 'k-')
axN.set_xlabel('Time ($t$)')
axN.set_ylabel('$N_{ph}$')
axN.set_title('Number of Phonons')
figN.savefig(dirpath + '/figures/quench_PhononNumber_aIBi:%.2f_P:%.2f.pdf' % (aIBi, P))

figPB, axPB = plt.subplots()
axPB.plot(tVec, PB_Vec, 'b-', tVec, P - PB_Vec, 'r--', tVec, nu(sParams[3]) + 0 * tVec, 'k-')
axPB.legend([r'$P_{B}$', r'$P_{I}$'])
axPB.set_xlabel('Time ($t$)')
axPB.set_ylabel('Momentum')
axPB.set_title('Phonon and Impurity Momentum')
figPB.savefig(dirpath + '/figures/quench_Phonon&ImpMomentum_aIBi:%.2f_P:%.2f.pdf' % (aIBi, P))

figPB, axPB = plt.subplots()
axPB.plot(tVec, P / (P - PB_Vec), 'k--')
axPB.set_xlabel('Time ($t$)')
axPB.set_ylabel('$M_{P}/M$')
axPB.set_title('Polaron Mass')
axPB.set_xlim([0, 1.5])
axPB.set_ylim([1, 10])
figPB.tight_layout()
figPB.savefig(dirpath + '/figures/quench_PolaronMass_aIBi:%.2f_P:%.2f.pdf' % (aIBi, P))


fig, axes = plt.subplots(nrows=1, ncols=2)
axes[0].plot(tVec, np.abs(S_Vec), 'k-')
axes[0].set_xlabel('Time ($t$)')
axes[0].set_ylabel(r'$\left|S(t)\right|$')
axes[0].set_title('Dynamical Overlap')


axes[1].plot(freqVec, A_Vec, 'k-')
axes[1].set_xlim([-100, 100])
# axes[1].set_ylim([0, 0.1])
axes[1].set_xlabel(r'Frequency ($\omega$)')
axes[1].set_ylabel(r'$A(\omega)$')
axes[1].set_title(r'Spectral Function')
fig.tight_layout()
fig.savefig(dirpath + '/figures/quench_DynOverlap&SpectFunction_aIBi:%.2f_P:%.2f.pdf' % (aIBi, P))

plt.show()
