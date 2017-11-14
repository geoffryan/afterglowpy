import sys
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import emcee as em
import grbpy as grb
import corner

lvars = [1,4,6,7,8,9]
labels_all = np.array([r"$\theta_{obs}$", r"$E_{iso}$", r"$\theta_j$",
                r"$\theta_w$", r"$n_0$", r"$p$", r"$\epsilon_e$",
                r"$\epsilon_B$", r"$\xi_N$", r"$d_L$"] )

def logpost(x, logprior, loglike, jetType, fluxArg, freeVars, bounds,
                tDat, nuDat, FnuDat, dFnuDat, opt=False):

    arg = fluxArg.copy()
    arg[freeVars] = x[:]
    arg[lvars] = np.power(10.0, arg[lvars])

    lp = logprior(jetType, bounds, *arg)

    if (x<bounds[:,0]).any() or (x>bounds[:,1]).any():
        lp = -np.inf

    if lp > -np.inf:
        lp += loglike(jetType, arg, tDat, nuDat, FnuDat, dFnuDat)

    print(str(x) + ": " + str(lp))

    if(opt):
        lp *= -1.0

    return lp

def chi2(jetType, arg, tDat, nuDat, FnuDat, dFnuDat):
    argtup = tuple(arg)
    Fnu = grb.fluxDensity(tDat, nuDat, jetType, *argtup)
    x = (Fnu-FnuDat) / dFnuDat
    x2 = (x*x).sum()

    return x2 

def logPriorFlat(jetType, bounds, *args):
    return 0.0


def logLikeChi2(jetType, arg, tDat, nuDat, FnuDat, dFnuDat):

    ch2 = chi2(jetType, arg, tDat, nuDat, FnuDat, dFnuDat)

    return -0.5*ch2

def plotChain(chain, fitPars):
    labels = labels_all[fitPars]
    ndim = len(fitPars)
    nwalkers = chain.shape[0]
    nsteps = chain.shape[1]

    for i in range(ndim):
        fig, ax = plt.subplots(1,1)
        for j in range(nwalkers):
            ax.plot(range(nsteps), chain[j,:,i], color='k', alpha=0.2)
        ax.set_xlabel("steps")
        ax.set_ylabel(labels[i])
        fig.tight_layout()
        fig.savefig("trace_" + str(fitPars[i]) + ".png")
        plt.close(fig)

    samples = chain[:,:,:].reshape((-1,ndim))

    fig = corner.corner(samples, labels=labels)
    fig.savefig("corner_all.png")
    plt.close(fig)

    if nwalkers > 20:
        for i in range(nsteps):
            fig = corner.corner(chain[:,i,:], labels=labels)
            fig.savefig("corner_{0:03d}.png".format(i))
            plt.close(fig)
    
def plot_curve(ax, t, FNU, Nt, Nnu, alpha=1.0):

    colors = ['k', 'b', 'g', 'r']

    for i in range(Nnu):
        ic = i % len(colors)
        ax.plot(t, Fnu[i*Nt:(i+1)*Nt], color=colors[ic], ls='-', 
                    marker='', alpha=alpha)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")


def plot_data(ax, t, FNU, FERR, Nt, Nnu):

    colors = ['k', 'b', 'g', 'r']

    for i in range(Nnu):
        ic = i % len(colors)
        ax.errorbar(t, FNU[i*Nt:(i+1)*Nt], FERR[i*Nt:(i+1)*Nt],
                    marker='', ls='', color=colors[ic])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")


if __name__ == "__main__":

    t0 = 1.0e5
    t1 = 1.0e7
    Nt = 30

    nu0 = 1.0e9
    nu1 = 1.0e18
    Nnu = 3

    t = np.logspace(np.log10(t0), np.log10(t1), num=Nt, base=10.0)
    nu = np.logspace(np.log10(nu0), np.log10(nu1), num=Nnu, 
                        base=10.0)

    N = Nt*Nnu

    T = np.empty(N)
    NU = np.empty(N)
    for i in range(Nnu):
        T[i*Nt:(i+1)*Nt] = t[:]
        NU[i*Nt:(i+1)*Nt] = nu[i]

    jetType = -1
    theta_obs = 0.3
    E_iso = 1.0e52
    theta_h = 0.1
    theta_h_wing = 0.2
    n0 = 1.0
    p = 2.5
    epsE = 0.1
    epsB = 0.01
    xiN = 1.0
    d_L = 1.0e28

    defaults = np.array([0.06, 1.0e53, 0.1, 0.2, 1.0, 2.2, 
                            0.1, 0.01, 1.0, 1.0e28])
    defaultPars = defaults.copy()
    defaultPars[lvars] = np.log10(defaultPars[lvars])

    errFac = 0.2

    #fitPars = [0,1,2,5,7]
    fitPars = [0,2,5]

    bounds = np.array([
                    [0.0, 0.5*np.pi],
                    [45.0, 57.0],
                    [0.01, 0.5*np.pi],
                    [0.0, 0.5*np.pi],
                    [-10.0, 10.0],
                    [1.0, 5.0],
                    [-10.0, 0.0],
                    [-10.0, 0.0],
                    [-10.0, 0.0],
                    [20, 40]])

    argtup = (theta_obs, E_iso, theta_h, theta_h_wing, n0, p, epsE,
            epsB, xiN, d_L)
    args = np.array(argtup)
    args[lvars] = np.log10(args[lvars])

    FNU = grb.fluxDensity(T, NU, jetType, *argtup)
    FERR = errFac * FNU

    ndim = len(fitPars)
    fitArgs = np.array(args)

    print("Truth: " + str(fitArgs[fitPars]))
    fitArgs[fitPars] = defaultPars[fitPars]

    #fitArgs[fitPars] *= 1.0 + errFac*(2*np.random.rand(ndim)-1)
    FNU *= 1.0 + errFac*np.random.randn(N)
   
    lpargs=(logPriorFlat, logLikeChi2, jetType, 
                fitArgs, fitPars, bounds[fitPars], T, NU, FNU, FERR, 
                False)
    """
    res = opt.minimize(logpost, fitArgs[fitPars],
                    args=lpargs,
                    bounds=bounds[fitPars])

    print(res)

    fitArgs[fitPars] = res.x[:]
    """

    nwalkers = 50
    nsteps = 100
    x0 = fitArgs[fitPars]
    noiseFac = 0.1
    p0 = [x0*(1+noiseFac*np.random.randn(ndim))
                                for i in range(nwalkers)]
    sampler = em.EnsembleSampler(nwalkers, ndim, logpost,
                                    args=lpargs)
    for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
        sys.stdout.write("\r{0:5.1%}    ".format(float(i)/nsteps))
    sys.stdout.write("\r{0:5.1%}    ".format(1.0))


    print("best: " + str(sampler.flatchain[i]))

    print("Plotting chain")
    plotChain(sampler.chain, fitPars)

    print("Plotting final ensemble")
    fig, ax = plt.subplots(1,1)
    for i in range(nwalkers):
        fluxArgs = fitArgs.copy()
        fluxArgs[fitPars] = sampler.chain[i,-1]
        fluxArgs[lvars] = np.power(10.0, fluxArgs[lvars])
        Fnu = grb.fluxDensity(T, NU, jetType, *fluxArgs)
        plot_curve(ax, t, Fnu, Nt, Nnu, alpha=0.2)
    plot_data(ax, t, FNU, FERR, Nt, Nnu)
    fig.savefig("lc_dist.png")
    plt.close(fig)

    print("Calculating best")
    i = np.argmax(sampler.flatlnprobability)
    fitArgs[fitPars] = sampler.flatchain[i]

    #fitArgs[fitPars] = args[fitPars]

    fluxArgs = fitArgs.copy()
    fluxArgs[lvars] = np.power(10.0, fluxArgs[lvars])

    Fnu = grb.fluxDensity(T, NU, jetType, *fluxArgs)

    print("Plotting Best")
    fig, ax = plt.subplots(1,1)

    plot_curve(ax, t, Fnu, Nt, Nnu)
    plot_data(ax, t, FNU, FERR, Nt, Nnu)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$F_\nu$ (mJy)")

    fig.savefig("lc_best.png")


