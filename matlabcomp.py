# from polrabi.basic import *
from polrabi.staticfm import *
# from polrabi.wf1 import *

gBB = 0.05
mI = 10
mB = 1
n0 = 1


# # basic.py test
# print (nu(gBB))
# print(ur(mI, mB))
# print(eB(20, mB))
# print(w(20, gBB, mB, n0))
# print((Wk(20, gBB, mB, n0)**2))


# staticfinmom.py test
w_rot = 3
P = 0.1
PB_v = 0.5
aIBi = -10
aSi_v = 0.3
k = 20
DP_v = 0.5

# print(Eup(P, PB_v, aIBi, aSi_v, mI, mB, n0))
# print(Edown(0, 0, aIBi, aSi_v, w_rot, mI, mB, n0))
# print(rMass(P, PB_v, mI))
# print(Bk(k, aIBi, aSi_v, gBB, mI, mB, n0))
# print(aSi0(gBB, mI, mB, n0))
# print(aSi(DP_v, gBB, mI, mB, n0))
# print(PB(DP_v, aIBi, gBB, mI, mB, n0))
# print(DP(DP_v, P, aIBi, gBB, mI, mB, n0))
# print(PMax(aIBi, gBB, mI, mB, n0))


# # wf1.py test
# Omega = 0.1
# w_rot = 3
# aIBi = -10
# aIBi_u = -10
# aIBi_d = 7
# aSi_v = 0.3

# print(cs_overlap(aIBi_u, aIBi_d, gBB, mI, mB, n0))
# print(qp_residue(aIBi, aSi_v, gBB, mI, mB, n0))
# print(Hspin(aIBi_u, aIBi_d, Omega, w_rot, gBB, mI, mB, n0))
