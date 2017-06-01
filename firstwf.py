from polrabi.wf1 import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import trapz, simps, cumtrapz


# # Initialization


matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

gBB = 0.05
mI = 10
mB = 1
n0 = 1


# # Rabi oscillation


# Omega = 0.1
# w_rot = 0
# aIBi_up = 30
# aIBi_down = 100
# a0 = np.array([[0.0], [1.0]], dtype=complex)

# tVals = np.linspace(0, 100, 1e3)
# a = np.zeros((2, tVals.size), dtype=complex)


# for idt, t in enumerate(tVals):
#     temp = wfcoeff(t, a0, aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0)
#     a[[0, 1], idt] = temp.reshape((2,))


# figR, axR = plt.subplots()
# axR.plot(tVals, abs(a[0, :])**2, 'k-')
# # ax.set_ylim([-50, 50])
# axR.set_xlabel('Time ($t$)')
# axR.set_ylabel(r'$P_{\uparrow}$')
# axR.set_title(r'Probability of being in $\uparrow$ state')
# # plt.show()


# # Atom transfer peak

Omega = 0.1
ts = np.pi / (2 * Omega)
aIBi_up = 100
aIBi_down = 100
a0 = np.array([[0.0], [1.0]], dtype=complex)

wVals = np.linspace(-15, 15, 1e3)
a = np.zeros((2, wVals.size), dtype=complex)

for idw, w in enumerate(wVals):
    temp = wfcoeff(ts, a0, aIBi_up, aIBi_down, Omega, w, gBB, mI, mB, n0)
    a[[0, 1], idw] = temp.reshape((2,))

figA, axA = plt.subplots()
axA.plot(wVals, np.abs(a[0, :])**2, 'k-')
# ax.set_ylim([-50, 50])
axA.set_xlabel(r'Pumping frequency ($\omega$)')
axA.set_ylabel(r'$P_{\uparrow}$')
axA.set_title(r'Probability of being in $\uparrow$ state')


p_up = np.abs(a[0, :])**2
p_down = np.abs(a[1, :])**2
mask = p_up > 0.1
print(trapz(p_up, x=wVals))
# print(p_up[mask])
# c = cumtrapz(p_up, x=wVals)
# print(c[mask[0:-1]])

# H = Hspin(aIBi_up, aIBi_down, Omega, w_rot, gBB, mI, mB, n0)
gam = 0.1
w21 = 0


def sakurai(w, t):
    return gam**2 / (gam**2 + (w - w21)**2 / 4) * np.sin(np.sqrt(gam**2 + (w - w21)**2 / 4) * t)**2


plt.show()
