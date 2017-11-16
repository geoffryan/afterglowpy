from __future__ import print_function
import sys
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import h5py as h5
import emcee as em
import corner
import grbpy as grb
import getOKCData as OKCdata

logVars = [1,4,6,7,8,9]
labelsAll = np.array([r"$\theta_{obs}$", r"$E_{iso}$", r"$\theta_j$",
                r"$\theta_w$", r"$n_0$", r"$p$", r"$\epsilon_e$",
                r"$\epsilon_B$", r"$\xi_N$", r"$d_L$"] )
day = 86400.0

bounds = np.array([[0.0, 0.5*np.pi], [45.0, 57.0], [0.01, 0.5*np.pi],
                    [0.0, 0.5*np.pi], [-10.0, 10.0], [1.0, 5.0], [-10.0, 0.0],
                    [-10.0, 0.0], [-10.0, 0.0], [20, 40]])

printLP = False

def logpost(x, logprior, loglike, jetType, fluxArg, fitVars, bounds,
                tDat, nuDat, FnuDat, dFnuDat, opt=False):

    X = fluxArg.copy()
    X[fitVars] = x[:]

    lp = logprior(x, jetType, X, fitVars, bounds)

    if lp > -np.inf:
        lp += loglike(jetType, X, tDat, nuDat, FnuDat, dFnuDat)

    if printLP:
        print(str(x) + ": " + str(lp))

    if(opt):
        lp *= -1.0

    return lp

def chi2(jetType, X, tDat, nuDat, FnuDat, dFnuDat):
    Y = getEvalForm(X)
    Fnu = grb.fluxDensity(tDat, nuDat, jetType, *Y)
    chi = (Fnu-FnuDat) / dFnuDat
    chi2 = (chi*chi).sum()

    return chi2 

def logPriorFlat(x, jetType, X, fitVars, bounds):

    lp = 0.0

    if bounds is not None:
        if (x<bounds[:,0]).any() or (x>bounds[:,1]).any():
            lp = -np.inf

    # Gaussian+Core, Wings must be larger than core.
    if jetType == 2 and X[3] < X[2]:
        lp = -np.inf

    return lp


def logLikeChi2(jetType, X, tDat, nuDat, FnuDat, dFnuDat):

    ch2 = chi2(jetType, X, tDat, nuDat, FnuDat, dFnuDat)

    return -0.5*ch2

def plotChain(chain, labels, fitVars):
    nwalkers = chain.shape[0]
    nsteps = chain.shape[1]
    ndim = chain.shape[2]

    for i in range(ndim):
        fig, ax = plt.subplots(1,1)
        for j in range(nwalkers):
            ax.plot(range(nsteps), chain[j,:,i], color='k', alpha=0.2)
        ax.set_xlabel("steps")
        ax.set_ylabel(labels[i])
        fig.tight_layout()
        fig.savefig("trace_" + str(fitVars[i]) + ".png")
        plt.close(fig)

    samples = chain[:,:,:].reshape((-1,ndim))

    fig = corner.corner(samples, labels=labels)
    fig.savefig("corner_all.png")
    plt.close(fig)

    #if nwalkers > 20:
    #    for i in range(nsteps):
    #        fig = corner.corner(chain[:,i,:], labels=labels)
    #        fig.savefig("corner_{0:03d}.png".format(i))
    #        plt.close(fig)
    
def plot_curve(ax, t, Fnu, color=None, alpha=1.0):

    colors = ['k', 'b', 'g', 'r']

    if color is None:
        color = 'k'

    ax.plot(t/day, Fnu, color=color, ls='-', marker='', alpha=alpha)


def plot_data(ax, t, Fnu, Ferr, inst):

    real = Fnu>0.0
    lim = Fnu<=0.0

    ax.errorbar(t[real]/day, Fnu[real], Ferr[real], color='b', ls='')
    ax.plot(t[lim]/day, 3*Ferr[lim], color='b', ls='', marker='v',
                                        mew=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (d)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")

def getEvalForm(X):
    Y = X.copy()
    Y[logVars] = np.power(10.0, Y[logVars])
    return Y

def getFitForm(Y):
    X = Y.copy()
    X[logVars] = np.log10(X[logVars])
    return X

def sample(X0, fitVars, jetType, bounds, data, nwalkers, nsteps, nburn, label):

    filename = label+".h5"
    ndim = len(fitVars)

    lpargs=(logPriorFlat, logLikeChi2, jetType, 
                X0, fitVars, bounds[fitVars],
                data[0], data[1], data[2], data[3], False)

    x0 = X0[fitVars]
    noiseFac = 0.02
    p0 = [x0*(1+noiseFac*np.random.randn(ndim))
                                for i in range(nwalkers)]

    nbuf = 10
    chainBuf = np.empty((nwalkers, nbuf, ndim))
    lnprobabilityBuf = np.empty((nwalkers, nbuf))


    f = h5.File(filename, "w")
    f.create_dataset("fitVars", data=fitVars)
    f.create_dataset("t", data=data[0])
    f.create_dataset("nu", data=data[1])
    f.create_dataset("Fnu", data=data[2])
    f.create_dataset("eFnu", data=data[3])
    f.create_dataset("inst", data=data[4].astype("S32"))
    f.create_dataset("X0", data=X0)
    f.create_dataset("jetType", data=np.array([jetType]))
    f.create_dataset("labels", data=labelsAll.astype("S32"))
    f.create_dataset("chain", (nwalkers, nsteps, ndim), dtype=np.float)
    f.create_dataset("lnprobability", (nwalkers, nsteps), dtype=np.float)
    f.create_dataset("steps_taken", data=np.array([1]))
    f.close()

    sampler = em.EnsembleSampler(nwalkers, ndim, logpost, args=lpargs)

    j=0
    k=0
    for i, result in enumerate(sampler.sample(p0, iterations=nsteps, 
                                                storechain=False)):
        chainBuf[:,j,:] = result[0]
        lnprobabilityBuf[:,j] = result[1]
        if j == nbuf-1:
            f = h5.File(filename, "a")
            f['chain'][:,k*nbuf:(k+1)*nbuf,:] = chainBuf[:,:,:]
            f['lnprobability'][:,k*nbuf:(k+1)*nbuf] = lnprobabilityBuf[:,:]
            f['steps_taken'][0] = (k+1)*nbuf
            f.close()
            k += 1
            j = 0
        else:
            j += 1
        sys.stdout.write("\r{0:5.1%}".format(float(i)/nsteps))
        sys.stdout.flush()
    f = h5.File(filename, "a")
    f['chain'][:,k*nbuf:k*nbuf+j,:] = chainBuf[:,:j,:]
    f['lnprobability'][:,k*nbuf:k*nbuf+j] = lnprobabiliyyBuf[:,:j]
    f.close()
    sys.stdout.write("\r{0:5.1%}\n".format(1.0))
    sys.stdout.flush()
    
    return sampler

def getDataTxt(datafile):

    t, nu, Fnu, Ferr, inst = np.loadtxt(datafile, unpack=True)
    return t, nu, Fnu, Ferr, inst

def getDataOKC(datafile):

    dat = OKCdata.OKCData("GW170817")
    tR, nuR, FnuR, eFnuR, instR = dat.getRadio()
    tX, nuX, FnuX, eFnuX, instX = dat.getXRay()

    realR = FnuR > 0.0
    tR = tR[realR]
    nuR = nuR[realR]
    FnuR = FnuR[realR]
    eFnuR = eFnuR[realR]
    instR = instR[realR]

    realX = FnuX > 0.0
    limX = FnuX <= 0.0
    bestLim = eFnuX[limX].min()
    keep = realX + (limX * (eFnuX<=1.1*bestLim))
    tX = tX[keep]
    nuX = nuX[keep]
    FnuX = FnuX[keep]
    eFnuX = eFnuX[keep]
    instX = instX[keep]

    tA = np.array([9. * 86400.0])
    nuA = np.array([nuX.mean()])
    FnuA = np.array([4.e-15 * 3.8e7])
    eFnuA = np.array([1.1e-15 * 3.8e7])
    instA = np.array(['Chandra'])

    t = np.concatenate((tR, tX, tA))
    nu = np.concatenate((nuR, nuX, nuA))
    Fnu = np.concatenate((FnuR, FnuX, FnuA))
    Ferr = np.concatenate((eFnuR, eFnuX, eFnuA))
    inst = np.concatenate((instR, instX, instA)).astype("S32")

    return t, nu, Fnu, Ferr, inst

def getPar(words, name):

    for line in words:
        if len(line) > 0 and line[0] != '#':
            if line[0] == name:
                return line[1]

    raise KeyError("Parameter ("+name+") not in parfile")

def getPars(words, name):

    for line in words:
        if len(line) > 0 and line[0] != '#':
            if line[0] == name:
                pars = []
                for val in line[1:]:
                    if val[0] != '#':
                        pars.append(val)
                    else:
                        break
                return pars

    raise KeyError("Parameter ("+name+") not in parfile")

def parseParfile(parfile):

    f = open(parfile, "r")
    lines = f.readlines()
    f.close()

    words = []
    for line in lines:
        words.append(line.split())

    label = getPar(words, "label")
    datafile = getPar(words, "datafile")
    nburn = int(getPar(words, "nburn"))
    nsteps = int(getPar(words, "nsteps"))
    nwalkers = int(getPar(words, "nwalkers"))
    fitVars = [int(x) for x in getPars(words, "fitVars")]
    jetType = int(getPar(words, "jetType"))
    thetaObs = float(getPar(words, "theta_obs"))
    Eiso = float(getPar(words, "E_iso_core"))
    thetaJ = float(getPar(words, "theta_h_core"))
    thetaW = float(getPar(words, "theta_h_wing"))
    n0 = float(getPar(words, "n_0"))
    p = float(getPar(words, "p"))
    epsE = float(getPar(words, "epsilon_E"))
    epsB = float(getPar(words, "epsilon_B"))
    xiN = float(getPar(words, "ksi_N"))
    dL = float(getPar(words, "d_L"))

    Y0 = np.array([thetaObs, Eiso, thetaJ, thetaW, n0, p, epsE, epsB, xiN, dL])
    X0 = getFitForm(Y0)

    return label, nwalkers, nburn, nsteps, fitVars, jetType, X0, datafile


def runFit(parfile):
    
    label, nwalkers, nburn, nsteps, fitVars, jetType, X0, datafile = parseParfile(parfile)

    dataExt = datafile.split(".")[-1]
    if dataExt == "json":
        T, NU, FNU, FERR, INST = getDataOKC(datafile)
    else:
        T, NU, FNU, FERR, INST = getDataTxt(datafile)

    data = (T, NU, FNU, FERR, INST)
    N = len(T)
    
    sampler = sample(X0, fitVars, jetType, bounds, data, nwalkers, nsteps,
                        nburn, label)

    print("Plotting chain")
    plotChain(sampler.chain, labelsAll[fitVars], fitVars)

    print("Plotting final ensemble")

    t0 = 1.0*day
    t1 = 300*day

    t = np.logspace(np.log10(t0), np.log10(t1), base=10.0, num=200)
    nu = np.empty(t.shape)
    nus = [6.0e9, 1.0e18]

    fig, ax = plt.subplots(1,1)

    X = X0.copy()
    for i in range(nwalkers):
        X[fitVars] = sampler.chain[i,-1,:]
        Y = getEvalForm(X)
        for v in nus:
            nu[:] = v
            Fnu = grb.fluxDensity(t, nu, jetType, *Y)
            plot_curve(ax, t, Fnu, alpha=0.2)
    plot_data(ax, T, FNU, FERR, INST)
    ax.set_xlim(t0/day, t1/day)
    ax.set_ylim(1.0e-9, 1.0e-1)
    fig.savefig("lc_dist.png")
    plt.close(fig)

    print("Calculating best")

    imax = np.argmax(sampler.flatlnprobability)
    print("best: " + str(sampler.flatchain[imax]))

    X1 = X0.copy()
    X1[fitVars] = sampler.flatchain[imax]
    Y1 = getEvalForm(X1)

    print("Plotting Best")

    fig, ax = plt.subplots(1,1)

    for v in nus:
        nu[:] = v
        Fnu = grb.fluxDensity(t, nu, jetType, *Y1)
        plot_curve(ax, t, Fnu)
    plot_data(ax, T, FNU, FERR, INST)
    ax.set_xlim(t0/day, t1/day)
    ax.set_ylim(1.0e-9, 1.0e-1)

    fig.savefig("lc_best.png")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please provide at least one parameter file")

    files = sys.argv[1:]
    for file in files:
        runFit(file)

