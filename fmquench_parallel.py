from polrabi.quench import *
from timeit import default_timer as timer
import os
from polrabi.staticfm import PCrit
import multiprocessing as mp
import itertools as it
from joblib import Parallel, delayed

# # Initialization


def dynamics(cParams, gParams, sParams):
    # takes parameters, performs dynamics, and outputs desired observables
    [P, aIBi] = cParams
    [kcutoff, dk, Ntheta, dtheta, tMax, dt] = gParams
    [mI, mB, n0, gBB] = sParams

    kVec = np.arange(dk, kcutoff, dk)
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

    # return data
    # dirpath = os.path.dirname(os.path.realpath(__file__))
    # np.save(dirpath + '/pdata/fmquench_aIBi:%.2f_P:%.2f.npy' % (aIBi, P), data)
    return data


if __name__ == "__main__":

    # set sParams

    mI = 1
    mB = 1
    n0 = 1
    gBB = (4 * np.pi / mB) * 0.05

    sParams = [mI, mB, n0, gBB]

    # set gParams

    kcutoff = 20
    dk = 0.05

    Ntheta = 10
    dtheta = np.pi / (Ntheta - 1)

    tMax = 4
    dt = 1e-5

    gParams = [kcutoff, dk, Ntheta, dtheta, tMax, dt]

    # create range of cParam values (P,aIBi)

    aIBi = -2
    Pc = PCrit(aIBi, gBB, mI, mB, n0)

    NPVals = 8
    PVals = np.linspace(0, 0.95 * Pc, NPVals)

    cParams_List = [[P, aIBi] for P in PVals]

    # create iterable over all tuples of function arguments for dynamics()

    paramsIter = zip(cParams_List, it.repeat(gParams), it.repeat(sParams))

    # compute data (parallel) - multiprocessing

    # start = timer()

    # with mp.Pool() as pool:
    #     # pool = mp.Pool()
    #     pool.starmap(dynamics, paramsIter)
    #     # pool.close()
    #     # pool.join()

    # end = timer()
    # print(end - start)

    # compute data (parallel) - joblib

    start = timer()

    num_cores = min(mp.cpu_count(), NPVals)
    print("Running on %d cores" % num_cores)
    results1 = Parallel(n_jobs=num_cores)(delayed(dynamics)(*p) for p in paramsIter)

    end = timer()
    print(end - start)

    # compute data (serial) - for loop

    # start = timer()

    # for z in paramsIter:
    #     dynamics(*z)

    # end = timer()
    # print(end - start)

    # compute data (serial) - starmap

    # start = timer()

    # for i in it.starmap(dynamics, paramsIter):
    #     i

    # end = timer()
    # print(end - start)
