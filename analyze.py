import sys
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import scipy.optimize as opt
import emcee
import corner
import grbpy as grb
import fit
import getOKCData as data

blue = (31.0/255, 119.0/255, 180.0/255)
orange = (255.0/255, 127.0/255, 14.0/255)
green = (44.0/255, 160.0/255, 44.0/255)
red = (214.0/255, 39.0/255, 40.0/255)
purple = (148.0/255, 103.0/255, 189.0/255)

fluxFunc = None

def plotTraces(chain, steps_taken, fitVars, labels, name):

    ndim = chain.shape[2]
    nwalkers = chain.shape[0]

    t = np.arange(steps_taken)

    for i in range(ndim):
        fig, ax = plt.subplots(1,1, figsize=fit.figsize)
        for j in range(nwalkers):
            trace = chain[j,:steps_taken,i]
            ax.plot(t, trace, alpha=0.1, color='k', ls='-')
        ax.set_xlabel("steps", fontsize=fit.labelsize)
        ax.set_ylabel(labels[fitVars[i]], fontsize=fit.labelsize)
        ax.tick_params(labelsize=fit.ticksize)
        fig.savefig("{0:s}_{1:d}.png".format(name, fitVars[i]))
        plt.close(fig)


def calcFluxStats(samples, fitVars, jetType, X0, t, nu):

    N = samples.shape[0]
    Fnu = np.empty((N,))
    T = np.array([t])
    NU = np.array([nu])

    X = X0.copy()

    for i in range(N):
        sys.stdout.write("\rCalculating flux {0:d} of {1:d}".format(i+1,N))
        sys.stdout.flush()
        X[fitVars] = samples[i]
        Y = fit.getEvalForm(jetType, X)
        Fnu[i] = fluxFunc(T, NU, jetType, *Y)
    sys.stdout.write("\n")

    mean = Fnu.mean()
    var = (Fnu*Fnu).mean()
    sig = np.sqrt(var - mean*mean)

    return mean, sig

def calcStats(chain):

    ndim = chain.shape[2]
    flatchain = chain.reshape((-1,ndim))

    tmeans = chain.mean(axis=0)
    tmom2 = (chain*chain).mean(axis=0)
    tvars = tmom2-tmeans*tmeans
    tquantiles = np.percentile(chain, (2.5,16,50,84,97.5), axis=0)
    
    means = tmeans.mean(axis=0)
    mom2 = tmom2.mean(axis=0)
    vars = mom2-means*means
    quantiles = np.percentile(flatchain, (2.5,16,50,84,97.5), axis=0)

    return means, vars, quantiles, tmeans, tvars, tquantiles

def calcAutocorr(chain, c=10):

    # Follow the suggestion of Daniel Foreman-Mackey, compute tau for
    # each walker then average. (http://dfm.io/posts/autocorr/)
    # Uses default parameters from emcee v 2.2.1

    #This is a straight copy of emcee's integration_time function,
    #but first taking an average of the walkers.

    size = 0.5*chain.shape[1]

    f = np.zeros((chain.shape[1], chain.shape[2]))
    for trace in chain:
        f += emcee.autocorr.function(trace, axis=0, fast=False)
    f /= chain.shape[0]
    m = [slice(None),] * len(f.shape)

    lo = 10
    hi = int(size / c)
    for M in np.arange(lo, hi).astype(int):
        m[0] = slice(1,M)
        taus = 1 + 2*np.sum(f[m], axis=0)
        if (taus > 1.0).all() and M > c * taus.max():
            return taus
    
    taus = np.empty(chain.shape[2])
    taus[:] = np.inf

    return taus

def plotStatsTimeSeries(tmeans, tsigs, tquantiles, fitVars, labels, name):

    steps = range(tmeans.shape[0])
    ndim = tmeans.shape[1]
    nq = tquantiles.shape[0]

    for i in range(ndim):
        fig, ax = plt.subplots(1,1, figsize=fit.figsize)
        ax.plot(steps, tmeans[:,i], color=blue, lw=2)
        ax.fill_between(steps, tmeans[:,i]-tsigs[:,i], tmeans[:,i]+tsigs[:,i],
                            color=blue, alpha=0.5)
        for j in range(nq):
            ax.plot(steps, tquantiles[j,:,i], color='k', ls='--')
        ax.set_xlabel("steps", fontsize=fit.labelsize)
        ax.set_ylabel(labels[i], fontsize=fit.labelsize)
        fig.savefig("{0:s}_{1:d}.png".format(name, fitVars[i]))
        plt.close(fig)

if __name__ == "__main__":
    

    file = sys.argv[1]

    print("Loading file.")

    f = h5.File(file, "r")
    steps_taken = f['steps_taken'][0]
    chain = f['chain'][...][:,:steps_taken,:]
    lnprobability = f['lnprobability'][...][:,:steps_taken]
    X0 = f['X0'][...]
    jetType = f['jetType'][0]
    fitVars = f['fitVars'][...]
    labels = f['labels'][...]
    T = f['t'][...]
    NU = f['nu'][...]
    FNU = f['Fnu'][...]
    FERR = f['eFnu'][...]
    UL = f['ul'][...]
    INST = f['inst'][...]
    f.close()
    
    if jetType == 3:
        bounds = fit.boundsCocoon
        fluxFunc = grb.fluxDensityCocoon
    else:
        bounds = fit.boundsJet
        fluxFunc = grb.fluxDensity

    #nburn = 14000
    #nburn = 100
    nburn = 1000
    nwalkers = chain.shape[0]
    nsteps = chain.shape[1]
    ndim = chain.shape[2]

    chainBurned = chain[:,nburn:,:]
    lnprobabilityBurned = lnprobability[:,nburn:]
    flatlnprobability = lnprobability.reshape((-1,))
    flatchain = chain.reshape((-1,ndim))
    flatlnprobabilityBurned = lnprobabilityBurned.reshape((-1,))
    flatchainBurned = chainBurned.reshape((-1,ndim))

    print("Calculating summary statistics")

    means, vars, quantiles, tmeans, tvars, tquantiles = calcStats(chain)
    tsigs = np.sqrt(tvars)
    print("Plotting summary statistic time series")
    plotStatsTimeSeries(tmeans, tsigs, tquantiles, fitVars, labels,
                            "trace_dist")

    print("Calculating burned-in summary statistics")

    means, vars, quantiles, tmeans, tvars, tquantiles = calcStats(chainBurned)
    sigs = np.sqrt(vars)
    tsigs = np.sqrt(tvars)
    print(means)
    print(vars)
    print(quantiles)

    print("Plotting burned-in summary statistic time series")
    plotStatsTimeSeries(tmeans, tsigs, tquantiles, fitVars, labels,
                            "trace_dist_noburn")
    
    print("Plotting walker traces")
    plotTraces(chain, steps_taken, fitVars, labels, "trace")
    print("Plotting burned-in walker traces")
    plotTraces(chainBurned, steps_taken-nburn, fitVars, labels, "trace_noburn")

    print("Autocorrelations:")
    taus = calcAutocorr(chain[:,:,:], 1)
    print(taus)
    taus = calcAutocorr(chain[:,:,:], 2)
    print(taus)
    taus = calcAutocorr(chain[:,:,:], 3)
    print(taus)
    taus = calcAutocorr(chain[:,:,:], 4)
    print(taus)
    taus = calcAutocorr(chain[:,:,:], 5)
    print(taus)
    taus = calcAutocorr(chain[:,nburn:,:], 1)
    print(taus)
    taus = calcAutocorr(chain[:,nburn:,:], 2)
    print(taus)
    taus = calcAutocorr(chain[:,nburn:,:], 3)
    print(taus)
    taus = calcAutocorr(chain[:,nburn:,:], 4)
    print(taus)
    taus = calcAutocorr(chain[:,nburn:,:], 5)
    print(taus)

    imax = np.argmax(flatlnprobability)
    x1 = flatchain[imax]
    print(x1)

    fig = corner.corner(flatchainBurned, labels=labels[fitVars],
                        quantiles=[0.16,0.50,0.84], truths=x1,
                        show_titles=True)
    fig.savefig("corner_noburn.png")
    plt.close(fig)

    fig = corner.corner(flatchain, labels=labels[fitVars],
                        quantiles=[0.16,0.50,0.84], truths=x1,
                        show_titles=True)
    fig.savefig("corner_all.png")
    plt.close(fig)

    N = len(T)

    X1 = X0.copy()
    X1[fitVars] = x1[:]

    lpargs=(fit.logPriorFlat, fit.logLikeChi2, jetType, 
                X1, fitVars, None, T, NU, FNU, FERR, 
                False)
    lpargsOpt=(fit.logPriorFlat, fit.logLikeChi2, jetType, 
                X1, fitVars, None, T, NU, FNU, FERR, 
                True)

    #f = h5.File(file, "r")
    #if "minChi2_val" in f.keys() and "minChi2_x" in f.keys():
    #    x2 = f['minChi2_x'][...]
    #    f.close()
    #else:
    #    f.close()
    #    fit.printLP = True
    #    res = opt.minimize(fit.logpost, x1, args=lpargsOpt, bounds=bounds[fitVars])
    #    fit.printLP = False
    #    x2 = res.x[:]

    #    f = h5.File(file, "a")
    #    f.create_dataset("minChi2_val", data=res.fun)
    #    f.create_dataset("minChi2_x", data=res.x)
    #    f.close()

    lp1 = flatlnprobability[imax]
    redChi21 = -2*lp1 / (ndim-1.0)
    print("best MCMC chi2/dof = " + str(redChi21))
    
    t0 = 1.0e0*fit.day
    t1 = 1.0e3*fit.day

    t = np.logspace(np.log10(t0), np.log10(t1), num=100, base=10.0)
    nu = np.empty(t.shape)
    FnuR = np.empty((nwalkers,t.shape[0]))
    FnuX = np.empty((nwalkers,t.shape[0]))

    nuX = NU[NU>1.0e16].mean()
    nuR = 6.0e9

    X = X0.copy()

    for i in range(nwalkers):
        sys.stdout.write("\rCalculating Fnu for walker " + str(i))
        sys.stdout.flush()
        X[fitVars] = chain[i,-1,:]
        Y = fit.getEvalForm(jetType, X)
        nu[:] = nuR
        FnuR[i,:] = fluxFunc(t, nu, jetType, *Y)
        nu[:] = nuX
        FnuX[i,:] = fluxFunc(t, nu, jetType, *Y)
    sys.stdout.write("\n")
    
    Y1 = fit.getEvalForm(jetType, X1)
    nu[:] = 6.0e9
    FnuR1 = fluxFunc(t, nu, jetType, *Y1)
    nu[:] = nuX.mean()
    FnuX1 = fluxFunc(t, nu, jetType, *Y1)

    #Y2 = fit.getEvalForm(X2)
    #nu[:] = 6.0e9
    #FnuR2 = grb.fluxDensity(t, nu, jetType, *Y2)
    #nu[:] = nuX.mean()
    #FnuX2 = grb.fluxDensity(t, nu, jetType, *Y2)
    
    fig, ax = plt.subplots(1,1, figsize=fit.figsize)

    for i in range(nwalkers):
        ax.plot(t/fit.day, FnuR[i], color='k', ls='-', alpha=0.1)
        ax.plot(t/fit.day, FnuX[i], color='k', ls='-', alpha=0.1)
    
    ax.plot(t/fit.day, FnuR1, color='k', ls='-', lw=2)
    ax.plot(t/fit.day, FnuX1, color='k', ls='-', lw=2)
    #ax.plot(t/fit.day, FnuR2, color='r', ls='-', lw=2)
    #ax.plot(t/fit.day, FnuX2, color='r', ls='-', lw=2)

    ax.set_xlim(t0/fit.day, t1/fit.day)
    ax.set_ylim(1.0e-9, 1.0e0)

    fit.plot_data(ax, T, NU, FNU, FERR, UL, INST)
    fig.savefig("lc_all.png")

    t_obs = 110.0 * fit.day

    ns = chainBurned.shape[1]
    step = 50
    samples = chainBurned[:,nsteps/2::step,:].reshape((-1,ndim))
    nsamps = samples.shape[0]

    mean110X, sig110X = calcFluxStats(samples, fitVars, jetType, X0,
                                        t_obs, nuX)
    mean110R, sig110R = calcFluxStats(samples, fitVars, jetType, X0,
                                        t_obs, nuR)

    print("110 day Radio Flux = {0:.3g} +/- {1:.3g} mJy (MCMC {2:.1g})".format(
                mean110R, sig110R, mean110R/np.sqrt(nsamps)))
    print("110 day X-Ray Flux = {0:.3g} +/- {1:.3g} mJy (MCMC {2:.1g})".format(
                mean110X, sig110X, mean110X/np.sqrt(nsamps)))

    sys.exit()

