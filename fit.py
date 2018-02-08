from __future__ import print_function
import os
import sys
import numpy as np
import scipy.optimize as opt
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py as h5
import emcee as em
import corner
import grbpy as grb

logVarsJet = [1,4,6,7,8,9]
labelsJetAll = np.array([r"$\theta_{obs}$", r"$E_{iso}$", r"$\theta_j$",
                r"$\theta_w$", r"$n_0$", r"$p$", r"$\epsilon_e$",
                r"$\epsilon_B$", r"$\xi_N$", r"$d_L$"] )

logVarsCocoon = [0,1,2,4,5,7,8,9,10]
labelsCocoonAll = np.array([r"$u_{max}$", r"$u_{min}$", r"$E_{inj}$",
                r"$k$", r"$M_{ej}$", r"$n_0$", r"$p$", r"$\epsilon_e$",
                r"$\epsilon_B$", r"$\xi_N$", r"$d_L$"] )
day = 86400.0

theta_min = 0.01

boundsJet = np.array([[0.0, 0.8], [45.0, 57.0], [theta_min, 0.5*np.pi],
                    [theta_min, 0.5*np.pi], [-10.0, 10.0], [2.0, 5.0],
                    [-4.0, 0.0], [-4.0, 0.0], [-4.0, 0.0], [20, 40]])
boundsCocoon = np.array([[-5, 3], [-5, 3], [45.0, 57.0], [0,8],
                    [-10, 0], [-10.0, 10.0], [2.0, 5.0], [-4.0, 0.0],
                    [-4.0, 0.0], [-4.0,0.0], [20.0, 40.0]])

printLP = False

figsize = (8,6)
labelsize = 24
legendsize = 18
ticksize = 18

blue = (31.0/255, 119.0/255, 180.0/255)
orange = (255.0/255, 127.0/255, 14.0/255)
green = (44.0/255, 160.0/255, 44.0/255)
red = (214.0/255, 39.0/255, 40.0/255)
purple = (148.0/255, 103.0/255, 189.0/255)
brown = (140.0/255, 86.0/255, 75.0/255)
pink = (227.0/255, 119.0/255, 194.0/255)
grey = (127.0/255, 127.0/255, 127.0/255)
yellow = (188.0/255, 189.0/255, 34.0/255)
teal = (23.0/255, 190.0/255, 207.0/255)

lightblue = (174.0/255, 199.0/255, 232.0/255)
lightorange = (255.0/255, 187.0/255, 120.0/255)
lightgreen = (152.0/255, 223.0/255, 138.0/255)
lightred = (255.0/255, 152.0/255, 150.0/255)
lightpurple = (197.0/255, 176.0/255, 213.0/255)
lightbrown = (196.0/255, 156.0/255, 148.0/255)
lightpink = (247.0/255, 182.0/255, 210.0/255)
lightgrey = (199.0/255, 199.0/255, 199.0/255)
lightyellow = (219.0/255, 219.0/255, 141.0/255)
lightteal = (23.0/255, 190.0/255, 207.0/255)

colors = [blue, orange, green, red, purple, brown, pink, grey, yellow, teal]
lightcolors = [lightblue, lightorange, lightgreen, lightred, lightpurple,
                    lightbrown, lightpink, lightgrey, lightyellow, lightteal]

def logpost(x, logprior, loglike, jetType, fluxArg, fitVars, bounds,
                tDat, nuDat, FnuDat, dFnuDat, opt=False):

    X = fluxArg.copy()
    X[fitVars] = x[:]

    lp = logprior(x, jetType, X, fitVars, bounds)

    if lp > -np.inf:
        #print(str(os.getpid()) + " - " + str(x))
        lp += loglike(jetType, X, tDat, nuDat, FnuDat, dFnuDat)

    if printLP:
        print(str(x) + ": " + str(lp))

    if(opt):
        lp *= -1.0

    return lp

def chi2(jetType, X, tDat, nuDat, FnuDat, dFnuDat):

    Y = getEvalForm(jetType, X)
    if jetType == 3:
        Fnu = grb.fluxDensityCocoon(tDat, nuDat, jetType, *Y)
    else:
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

    #Physics
    if jetType != 3 and np.isfinite(lp):
        lp += np.log(np.sin(X[0]))

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

    fig = corner.corner(samples, labels=labels, quantiles=[0.16,0.5,0.84],
                            show_titles=True)
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


def plot_data(ax, t, nu, Fnu, Ferr, ul, inst, spec=False, legend=True):

    cmapR = mpl.cm.get_cmap('Greens')
    cmapO = mpl.cm.get_cmap('Blues')
    cmapX = mpl.cm.get_cmap('Purples')

    insts = np.unique(inst)
    for instrument in insts:
        ind = inst==instrument
        nus = np.unique(nu[ind])
        for v in nus:
            ind2 = v==nu[ind]
            myt = t[ind][ind2]
            mynu = nu[ind][ind2]
            myFnu = Fnu[ind][ind2]
            myFerr = Ferr[ind][ind2]
            myul = ul[ind][ind2]

            if v < 1.0e11:
                nuRhi = 1.0e11
                nuRlo = 1.0e9
                mycolor=cmapR(np.log(v/nuRhi)/np.log(nuRlo/nuRhi))
                label = "{0:.1f} GHz".format(v/1.0e9)
            elif v < 1.0e15:
                nuOhi = 1.0e14
                nuOlo = 1.0e15
                mycolor=cmapO(np.log(v/nuOhi)/np.log(nuOlo/nuOhi))
                label = "i"
            elif v < 1.0e20:
                if instrument == 'NuSTAR':
                    mycolor=cmapX(0.2)
                    label = instrument
                elif instrument == 'Swift':
                    mycolor=cmapX(0.5)
                    label = instrument
                elif instrument == 'Chandra':
                    mycolor=cmapX(0.8)
                    label = instrument
                else:
                    mycolor=cmapX(1.0)
                    label = ''
            else:
                mycolor='k'
                label = ''
            real = myul <= 0.0
            lim = myul > 0.0
            if lim.any():
                if not spec:
                    ax.plot(myt[lim]/day, (myul*myFerr)[lim],
                            marker='v', color=mycolor, mew=0.0, ls='',
                            ms=10, label=label)
                else:
                    ax.plot(mynu[lim], (myul*myFerr)[lim],
                            marker='v', color=mycolor, mew=0.0, ls='',
                            ms=10, label=label)
            if real.any():
                if not spec:
                    ax.errorbar(myt[real]/day, myFnu[real], myFerr[real],
                            marker='', color=mycolor, ls='',
                            lw=2, label=label)
                else:
                    ax.errorbar(mynu[real], myFnu[real], myFerr[real],
                            marker='', color=mycolor, ls='',
                            lw=2, label=label)

    formatAxes(ax, legend, spec)


def formatAxes(ax, legend=True, spec=False, loc='upper left'):

    #if not spec:
    #    ax.axvline(110, lw=4, ls='--', color='grey')

    if legend:
        plt.legend(loc=loc, fontsize=legendsize)
    ax.set_xscale("log")
    ax.set_yscale("log")
    if not spec:
        ax.set_xlabel(r"$t$ (d)", fontsize=labelsize)
    else:
        ax.set_xlabel(r"$\nu$ (Hz)", fontsize=labelsize)
    ax.set_ylabel(r"$F_\nu$ (mJy)", fontsize=labelsize)
    ax.tick_params(labelsize=ticksize)
    if not spec:
        ax.set_xlim(1.0, 1.0e3)
    else:
        ax.set_xlim(1.0e6, 1.0e20)
    ax.set_ylim(1.0e-9, 1.0e0)
    ax.get_figure().tight_layout()

def dumpLCTxt(t, nu, Fnu, jetType, Y, name):

    f = open(name, "w")
    f.write("# nu={0:.6g} Hz\n".format(nu))
    f.write("# " + str(jetType) + " " + " ".join([str(y) for y in Y]) + "\n")
    f.write("# t(d) Fnu(mJy)\n")
    for i in range(t.shape[0]):
        f.write("{0:.8g} {1:.8g}\n".format(t[i]/day, Fnu[i]))
    f.close()


def getEvalForm(jetType, X):
    Y = X.copy()
    if jetType == 3:
        Y[logVarsCocoon] = np.power(10.0, Y[logVarsCocoon])
    else:
        Y[logVarsJet] = np.power(10.0, Y[logVarsJet])
    return Y

def getFitForm(jetType, Y):
    X = Y.copy()
    if jetType == 3:
        X[logVarsCocoon] = np.log10(X[logVarsCocoon])
    else:
        X[logVarsJet] = np.log10(X[logVarsJet])
    return X

def sample(X0, fitVars, jetType, bounds, data, nwalkers, nsteps, nburn, label,
            threads, restart=False):

    filename = label+".h5"
    ndim = len(fitVars)

    lpargs=(logPriorFlat, logLikeChi2, jetType, 
                X0, fitVars, bounds[fitVars],
                data[0], data[1], data[2], data[3], False)

    nbuf = 10
    chainBuf = np.empty((nwalkers, nbuf, ndim))
    lnprobabilityBuf = np.empty((nwalkers, nbuf))

    if jetType == 3:
        varLabels = labelsCocoonAll
    else:
        varLabels = labelsJetAll

    if not restart:
        f = h5.File(filename, "w")
        f.create_dataset("fitVars", data=fitVars)
        f.create_dataset("t", data=data[0])
        f.create_dataset("nu", data=data[1])
        f.create_dataset("Fnu", data=data[2])
        f.create_dataset("eFnu", data=data[3])
        f.create_dataset("ul", data=data[4])
        f.create_dataset("inst", data=data[5].astype("S32"))
        f.create_dataset("X0", data=X0)
        f.create_dataset("nburn", data=np.array([nburn]))
        f.create_dataset("jetType", data=np.array([jetType]))
        f.create_dataset("labels", data=varLabels.astype("S32"))
        f.create_dataset("chain", (nwalkers, nsteps, ndim), dtype=np.float)
        f.create_dataset("lnprobability", (nwalkers, nsteps), dtype=np.float)
        f.create_dataset("steps_taken", data=np.array([1]))
        f.create_dataset("threads", data=np.array([threads]))
        f.close()
        steps_taken = 0
    else:
        f = h5.File(filename, "r")
        steps_taken = f['steps_taken'][0]
        p0 = f['chain'][:,steps_taken-1,:]
        f.close()

    if steps_taken == 0:
        x0 = X0[fitVars]
        noiseFac = 1.0e-4
        p0 = np.array([x0*(1+noiseFac*np.random.randn(ndim))
                                    for i in range(nwalkers)])
        for i in range(ndim):
            if x0[i] == 0.0:
                p0[:,i] = noiseFac * np.random.randn(nwalkers)

    #for i in range(nwalkers):
    #    p = logpost(p0[i], *lpargs)
    #    print(str(i) + ": " + str(p) + " - " + str(p0[i]))

    sampler = em.EnsembleSampler(nwalkers, ndim, logpost, args=lpargs,
                                    threads=threads)

    j=0
    k=0
    for i, result in enumerate(sampler.sample(
                        p0, iterations=nsteps-steps_taken, storechain=False)):
        chainBuf[:,j,:] = result[0]
        lnprobabilityBuf[:,j] = result[1]
        if j == nbuf-1:
            k0 = steps_taken + k*nbuf
            f = h5.File(filename, "a")
            f['chain'][:,k0:k0+nbuf,:] = chainBuf[:,:,:]
            f['lnprobability'][:,k0:k0+nbuf] = lnprobabilityBuf[:,:]
            f['steps_taken'][0] = k0+nbuf
            f.close()
            k += 1
            j = 0
        else:
            j += 1
        sys.stdout.write("\r{0:5.1%}".format(float(i+steps_taken)/nsteps))
        sys.stdout.flush()
    k0 = steps_taken + k*nbuf
    f = h5.File(filename, "a")
    f['chain'][:,k0:k0+j,:] = chainBuf[:,:j,:]
    f['lnprobability'][:,k0:k0+j] = lnprobabilityBuf[:,:j]
    f['steps_taken'][0] = k0+j
    f.close()
    sys.stdout.write("\r{0:5.1%}\n".format(1.0))
    sys.stdout.write("Mean Acceptance Fraction: {0:.3f}\n".format(np.mean(
                        sampler.acceptance_fraction)))
    sys.stdout.flush()
    
    return sampler

def getDataTxt(datafile):

    t, nu, Fnu, Ferr, ul = np.loadtxt(datafile, unpack=True,
                                            usecols=[0,1,2,3,4])
    inst = np.loadtxt(datafile, unpack=True, usecols=[5], dtype='S32')
    return t, nu, Fnu, Ferr, ul, inst

def getDataOKC(datafile):

    import getOKCData as OKCdata
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

    ul = np.empty(t.shape)
    ul[Fnu>0.0] = 0.0
    ul[Fnu<=0.0] = 3.0

    return t, nu, Fnu, Ferr, ul, inst

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
    threads = int(getPar(words, "threads"))
    fitVars = [int(x) for x in getPars(words, "fitVars")]
    jetType = int(getPar(words, "jetType"))
    if jetType == 3:
        umax = float(getPar(words, "u_max"))
        umin = float(getPar(words, "u_min"))
        Ei = float(getPar(words, "E_i"))
        k = float(getPar(words, "k"))
        Mej = float(getPar(words, "M_ej"))
        n0 = float(getPar(words, "n_0"))
        p = float(getPar(words, "p"))
        epsE = float(getPar(words, "epsilon_E"))
        epsB = float(getPar(words, "epsilon_B"))
        ksiN = float(getPar(words, "ksi_N"))
        dL = float(getPar(words, "d_L"))

        Y0 = np.array([umax, umin, Ei, k, Mej, n0, p, epsE, epsB, ksiN, dL])
    else:
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
    X0 = getFitForm(jetType, Y0)

    return label, nwalkers, nburn, nsteps, fitVars, jetType, X0, threads, datafile


def runFit(parfile):

    parExt = parfile.split(".")[-1]
    if parExt == "h5":
        restart = True
        label = ".".join(parfile.split(".")[:-1])
        f = h5.File(parfile, "r")
        T = f['t'][...]
        NU = f['nu'][...]
        FNU = f['Fnu'][...]
        FERR = f['eFnu'][...]
        UL = f['ul'][...]
        INST = f['inst'][...]
        nwalkers = f['chain'].shape[0]
        nsteps = f['chain'].shape[1]
        ndim = f['chain'].shape[2]
        X0 = f['X0'][...]
        jetType = f['jetType'][0]
        fitVars = f['fitVars'][...]
        nburn = f['nburn'][0]
        threads = f['threads'][0]
        f.close()

    else:
        restart = False
        label, nwalkers, nburn, nsteps, fitVars, jetType, X0, threads, datafile = parseParfile(parfile)
        ndim = len(fitVars)

        dataExt = datafile.split(".")[-1]
        if dataExt == "json":
            T, NU, FNU, FERR, UL, INST = getDataOKC(datafile)
        else:
            T, NU, FNU, FERR, UL, INST = getDataTxt(datafile)

    data = (T, NU, FNU, FERR, UL, INST)
    N = len(T)

    if jetType == 3:
        bounds = boundsCocoon
        varLabels = labelsCocoonAll
    else:
        bounds = boundsJet
        varLabels = labelsJetAll
    
    sampler = sample(X0, fitVars, jetType, bounds, data, nwalkers, nsteps,
                        nburn, label, threads, restart=restart)

    print("Plotting chain")
    f = h5.File(label+".h5", "r")
    chain = f['chain'][...]
    lnprobability = f['lnprobability'][...]
    flatchain = chain.reshape((-1,ndim))
    flatlnprobability = lnprobability.reshape((-1,))
    f.close()
    plotChain(chain, varLabels[fitVars], fitVars)


    print("Plotting final ensemble")

    t0 = 1.0*day
    t1 = 300*day

    t = np.logspace(np.log10(t0), np.log10(t1), base=10.0, num=200)
    nu = np.empty(t.shape)
    nus = [6.0e9, 1.0e18]

    fig, ax = plt.subplots(1,1, figsize=figsize)

    if jetType == 3:
        flux = grb.fluxDensityCocoon
    else:
        flux = grb.fluxDensity

    X = X0.copy()
    for i in range(nwalkers):
        X[fitVars] = chain[i,-1,:]
        Y = getEvalForm(jetType, X)
        for v in nus:
            nu[:] = v
            Fnu = flux(t, nu, jetType, *Y)
            plot_curve(ax, t, Fnu, alpha=0.2)
    plot_data(ax, T, NU, FNU, FERR, UL, INST)
    fig.savefig("lc_dist.png")
    plt.close(fig)

    print("Calculating best")

    imax = np.argmax(flatlnprobability)
    print("best: " + str(flatchain[imax]))

    X1 = X0.copy()
    X1[fitVars] = flatchain[imax]
    Y1 = getEvalForm(jetType, X1)

    print("Plotting Best")

    fig, ax = plt.subplots(1,1, figsize=figsize)

    for v in nus:
        nu[:] = v
        Fnu = flux(t, nu, jetType, *Y1)
        plot_curve(ax, t, Fnu)
    plot_data(ax, T, NU, FNU, FERR, UL, INST)

    fig.savefig("lc_best.png")


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Please provide at least one parameter file")

    files = sys.argv[1:]
    for file in files:
        runFit(file)

