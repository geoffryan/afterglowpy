import numpy as np
import matplotlib.pyplot as plt
import grbpy as grb


def testLatRes(t, nu, jetType, specType, Y, rtol, lRs):

    Nres = len(lRs)
    Fnu = np.empty((Nres, t.shape[0]))

    print("latRes")
    for i,lR in enumerate(lRs):
        print(lR)
        Fnu[i,:] = grb.fluxDensity(t, nu, jetType, specType, *Y, 
                                    latRes=lR, rtol=rtol)
    FnuE = grb.fluxDensity(t, nu, jetType, specType, *Y, 
                                    latRes=2*lRs[-1], rtol=rtol)

    fig, ax = plt.subplots(2,1)
    for i in range(Nres):
        ax[0].plot(t, Fnu[i])
    ax[0].plot(t, FnuE, 'k')
    for i in range(Nres):
        ax[1].plot(t, np.fabs(Fnu[i]/FnuE-1.0))

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")

    print("saving latRes.png")
    fig.savefig("latRes.png", dpi=300)

def testRTol(t, nu, jetType, specType, Y, latRes, rtols):
    
    Nrt = len(rtols)

    Fnu = np.empty((Nrt, t.shape[0]))

    print("rtol")
    for i,rt in enumerate(rtols):
        print(rt)
        Fnu[i,:] = grb.fluxDensity(t, nu, jetType, specType, *Y, 
                                    latRes=latRes, rtol=rt)
    FnuE = grb.fluxDensity(t, nu, jetType, specType, *Y, 
                                    latRes=latRes, rtol=0.001*rtols[-1])

    fig, ax = plt.subplots(2,1)
    for i in range(Nrt):
        ax[0].plot(t, Fnu[i])
        ax[1].plot(t, np.fabs(Fnu[i]/FnuE-1.0))
    ax[0].plot(t, FnuE, 'k')

    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")

    print("saving rtol.png")
    fig.savefig("rtol.png", dpi=300)


if __name__ == "__main__":
    t = np.logspace(-1, 3, 100) * grb.day2sec
    nu = np.empty(t.shape)
    nu[:] = 1.0e18

    jetType = 0
    specType = 0
    Y = np.array([0.5, 1.0e53, 0.05, 0.4, 1.0e-3, 2.2, 0.1, 0.01, 1.0, 1.0e26])
    rtol = 1.0e-6

    Nres = 10
    lRs = range(1, Nres+1)
    testLatRes(t, nu, jetType, specType, Y, rtol, lRs)

    latRes = 5

    rtols = [1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8]
    
    testRTol(t, nu, jetType, specType, Y, latRes, rtols)
