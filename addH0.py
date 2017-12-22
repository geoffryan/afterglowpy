import sys
import math
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate
import emcee
import corner
import grbpy as grb
import fit

optimizeMAP = False
thin = 1

blue = (31.0/255, 119.0/255, 180.0/255)
orange = (255.0/255, 127.0/255, 14.0/255)
green = (44.0/255, 160.0/255, 44.0/255)
red = (214.0/255, 39.0/255, 40.0/255)
purple = (148.0/255, 103.0/255, 189.0/255)

cornerLabelsize = 24
cornerTitlesize = 18

fluxFunc = None

def weightPlanck(flatchain):

    thv = flatchain[:,0]
    cthv = np.cos(thv)

    #From Fitting Figure 3 of LIGO Standard Siren Paper
    cv0 = 0.985
    sig = 0.070

    x = (cthv-cv0)/sig
    w = np.exp(-0.5*x*x)

    return w

def weightSHoES(flatchain):

    thv = flatchain[:,0]
    cthv = np.cos(thv)

    cv0 = 0.909
    sig = 0.068

    x = (cthv-cv0)/sig
    w = np.exp(-0.5*x*x)

    return w

def weightLIGO(flatchain):

    thv = flatchain[:,:,0]
    cthv = np.cos(thv)

    cv0 = 0.958
    sig = 0.266

    x = (cthv-cv0)/sig
    w = np.exp(-0.5*x*x)

    return w

def weightEtot(flatchain, jetType, X0, fitVars):

    ndim = flatchain.shape[1]

    lEtot = f_Etot(flatchain, jetType, X0, fitVars)
    lEsolar = np.log10(grb.Msun * grb.c*grb.c)

    w = np.zeros((flatchain.shape[0],))
    w[lEtot < lEsolar] = 1.0

    return w

def f_Etot(flatchain, jetType, X0, fitVars):
    if jetType == 0:
        if 1 in fitVars:
            lE0 = flatchain[:,1]
        else:
            lE0 = np.atleast_1d(X0[1])
        if 2 in fitVars:
            thC = flatchain[:,2]
        else:
            thC = np.atleast_1d(X0[2])
        if 3 in fitVars:
            thW = flatchain[:,3]
        else:
            thW = np.atleast_1d(X0[3])
        
        x = thW/thC
        w = np.exp(-0.5*x*x)
        # thInt = int dphi dtheta sin(th) e^(-thC^2/2*thW^2) / 4*pi
        thInt = thC*thC * ((6-2*thC*thC)*(1-w) + w*thW*thW) / 12.0

        #Etot = 2 * Etot * thInt
        lEtot = lE0 + np.log10(2 * thInt)
    
    elif jetType == 3:
        if 0 in fitVars:
            lumax = flatchain[:,0]
        else:
            lumax = np.atleast_1d(X0[0])
        if 1 in fitVars:
            lumin = flatchain[:,1]
        else:
            lumin = np.atleast_1d(X0[1])
        if 2 in fitVars:
            lEi = flatchain[:,2]
        else:
            lEi = np.atleast_1d(X0[2])
        if 3 in fitVars:
            k = flatchain[:,3]
        else:
            k = np.atleast_1d(X0[3])
        if 4 in fitVars:
            lMej = flatchain[:,4]
        else:
            lMej = np.atleast_1d(X0[4])

        Einj = np.power(10.0, lEi - k*lumin) - np.power(10.0, lEi - k*lumax)
        Mej = grb.Msun * np.power(10.0, lMej) * grb.c*grb.c
        gmax = np.sqrt(1.0 + np.power(10.0, 2*lumax))

        # Etot = (gmax-1) * Mej + Einj
        lEtot = np.log10((gmax-1) * Mej + Einj)

    return lEtot

def f_thV_o_thC(flatchain, jetType, X0, fitVars):

    if 0 in fitVars:
        thV = flatchain[:,0]
    else:
        thV = np.atleast_1d(X0[0])

    if 2 in fitVars:
        thC = flatchain[:,2]
    else:
        thC = np.atleast_1d(X0[2])

    return thV / thC

def f_thW_o_thC(flatchain, jetType, X0, fitVars):

    if 3 in fitVars:
        thW = flatchain[:,3]
    else:
        thW = np.atleast_1d(X0[3])

    if 2 in fitVars:
        thC = flatchain[:,2]
    else:
        thC = np.atleast_1d(X0[2])

    return thW / thC


def dEjet_gaussian(theta, theta_j):
    return math.sin(theta)*math.exp(-0.5*theta*theta/(theta_j*theta_j))

def calcEjet_quick(jetType, x, X0, fitVars):
    if jetType == 0:
        theta_j = Y[2]
        theta_w = Y[3]
        Eiso_core = Y[1]
        res = integrate.quad(dEjet_gaussian, 0.0, theta_w, args=(theta_j,))
        thetaI = res[0]
        # Eiso/4pi * (integral over phi) * (integral over theta)
        # Eiso/4pi * 2pi * ...
        # 0.5*Eiso * (integral over theta)
        Ejet = 0.5 * Eiso_core * thetaI
    elif jetType == 1:
        Eiso = Y[1]
        theta_j = Y[2]
        theta_w = Y[3]

        # int_{theta_j}^{theta_w} (theta_j/theta)^2 dtheta
        thetaI = theta_j*(1.0 - theta_j/theta_w)
        Ejet = 0.5 * Eiso * (1-math.cos(theta_j) + theta_I)

    elif jetType == 3:
        umax = Y[0]
        umin = Y[1]
        Ei = Y[2]
        k = Y[3]
        Mi = Y[4] * grb.Msun * grb.c*grb.c

        gmax = math.sqrt(1+umax*umax)
        EiTot = Ei*(math.pow(umin,-k) - math.pow(umax,-k))
        Ejet = (gmax-1) * Mi + EiTot
    else:
        # Tophat (-1)
        Eiso = Y[1]
        theta_j = Y[2]
        Ejet = 0.5 * Eiso * (1-math.cos(theta_j))

    return Ejet



def analyze(flatchain, flatlp, weights, jetType, X0, fitVars, labels, data, 
                label, optimizeMAP=False, derivedVals=None):
    
    ndim = flatchain.shape[1]
    print("Analyzing {0:s}".format(label))
    #nwalkers = chain.shape[0]
    #nsteps = chain.shape[1]
    #ndim = chain.shape[2]
    #flatchain = chain.reshape((-1,ndim))
    #flatlp = lp.reshape((-1,))

    print("Calculating quantiles")
    for i in range(ndim):
        q = corner.quantile(flatchain[:,i], [0.025, 0.16, 0.5, 0.84, 0.975],
                                weights)
        print("{0:d}: {1:.5g} {2:.5g} {3:.5g} {4:.5g} {5:.5g}".format(
                i, q[0], q[1], q[2], q[3], q[4]))
        print("       68%: ({0:.2g} {1:.2g}) 95%: ({2:.2g} {3:.2g})".format(
                q[3]-q[2], q[1]-q[2], q[4]-q[2], q[0]-q[2]))

    print("Finding MAP")
    if weights is None:
        imax = np.argmax(flatlp)
    else:
        imax = np.argmax(flatlp + np.log(weights))
    xMAP0 = flatchain[imax]
    print("MAP0: ", xMAP0, flatlp[imax])

    if optimizeMAP:
        print("Optimizing MAP")
        lpargs=(fit.logPriorFlat, fit.logLikeChi2, jetType, X0, fitVars,
                None, data[0], data[1], data[2], data[3], True)
        if jetType == 3:
            bounds = fit.boundsCocoon[fitVars]
        else:
            bounds = fit.boundsJet[fitVars]

        fit.printLP = True
        res = opt.minimize(fit.logpost, xMAP0, args=lpargs, bounds=bounds)
        fit.printLP = False
        xMAP = res.x[:]
    else:
        xMAP = xMAP0
    print("MAP: ", xMAP)


    NDataTotal = len(data[0])
    NUL = len(data[0][data[4]==0.0])

    XMAP = X0.copy()
    XMAP[fitVars] = xMAP[:]

    T = data[0]
    NU = data[1]
    FNU = data[2]
    FERR = data[3]
    UL = data[4]
    N = len(T)
    
    obs = (UL==0.0)
    To = T[obs]
    NUo = NU[obs]
    FNUo = FNU[obs]
    FERRo = FERR[obs]
    No = len(To)

    try:
        chi2MAP = fit.chi2(jetType, XMAP, T, NU, FNU, FERR) 
        redchi2MAP = chi2MAP / (N - ndim)
        chi2MAPo = fit.chi2(jetType, XMAP, To, NUo, FNUo, FERRo) 
        redchi2MAPo = chi2MAPo / (No - ndim)
        print("MAP(obs+ul) chi2={0:.3g} chi2/dof={1:.3g}".format(
                                                        chi2MAP, redchi2MAP))
        print("MAP(obs)    chi2={0:.3g} chi2/dof={1:.3g}".format(
                                                        chi2MAPo, redchi2MAPo))
        t0 = 1.0e0*fit.day
        t1 = 1.0e3*fit.day

        t = np.logspace(np.log10(t0), np.log10(t1), num=100, base=10.0)
        nu = np.empty(t.shape)

        nuX = NU[NU>1.0e16].mean()
        nuR = 6.0e9

        Y = fit.getEvalForm(jetType, XMAP)
        nu[:] = nuR
        FnuRMAP = fluxFunc(t, nu, jetType, *Y)
        nu[:] = nuX
        FnuXMAP = fluxFunc(t, nu, jetType, *Y)

        fig, ax = plt.subplots(1,1, figsize=fit.figsize)
        ax.plot(t/fit.day, FnuRMAP, color='k', ls='-', lw=2)
        ax.plot(t/fit.day, FnuXMAP, color='k', ls='-', lw=2)
        fit.plot_data(ax, *data, legend=False)
        fig.savefig("{0:s}_lc_MAP.png".format(label))
        plt.close(fig)

        fit.dumpLCTxt(t, nuR, FnuRMAP, jetType, Y, 
                        "{0:s}_lc_MAP_radio.txt".format(label))
        fit.dumpLCTxt(t, nuX, FnuXMAP, jetType, Y, 
                        "{0:s}_lc_MAP_xray.txt".format(label))

        nu = np.logspace(3, 20, num=t.shape[0], base=10.0)
        t[:] = 17*fit.day
        Fnu17dMAP = fluxFunc(t, nu, jetType, *Y)
        t[:] = 110*fit.day
        Fnu110dMAP = fluxFunc(t, nu, jetType, *Y)

        fig, ax = plt.subplots(1,1, figsize=fit.figsize)
        ax.plot(nu, Fnu17dMAP, color='k', ls='-', lw=2)
        ax.plot(nu, Fnu110dMAP, color='k', ls='-', lw=2)
        fit.plot_data(ax, *data, legend=False, spec=True)
        fig.savefig("{0:s}_spec_MAP.png".format(label))
        plt.close(fig)

    except ValueError:
        print("Value Error in MAP calculation??")
    

    print("Making Corner Plot")
    fig = corner.corner(flatchain, labels=labels[fitVars],
                        quantiles=[0.16,0.50,0.84], truths=xMAP,
                        weights=weights, show_titles=True,
                        label_kwargs={'fontsize': cornerLabelsize},
                        title_kwargs={'fontsize': cornerTitlesize})
    fig.savefig("{0:s}_corner.png".format(label))
    plt.close(fig)

    if derivedVals is not None:
        print("Calculating Derived Vals")
        nval = len(derivedVals)
        labelsFull = np.empty((ndim+nval,), dtype=labels.dtype)
        flatchainFull = np.empty((flatchain.shape[0], nval+ndim),
                                    dtype=flatchain.dtype)
        truthFull = np.empty((ndim+nval,))
        labelsFull[:ndim] = labels[fitVars][:]
        flatchainFull[:,:ndim] = flatchain[:,:]
        truthFull[:ndim] = xMAP[:]
        truth = np.array([xMAP])
        for i in range(nval):
            labelsFull[ndim+i] = derivedVals[i]['label']
            f = derivedVals[i]['f']
            flatchainFull[:,ndim+i] = f(flatchain, jetType, X0, fitVars)
            truthFull[ndim+i] = f(truth, jetType, X0, fitVars)
        print("Making Full Corner Plot")
        fig = corner.corner(flatchainFull, labels=labelsFull,
                            quantiles=[0.16,0.50,0.84], truths=truthFull,
                            weights=weights, show_titles=True,
                            label_kwargs={'fontsize': cornerLabelsize},
                            title_kwargs={'fontsize': cornerTitlesize})
        fig.savefig("{0:s}_corner_full.png".format(label))
        plt.close(fig)



if __name__ == "__main__":
    

    file = sys.argv[1]

    print("Loading file.")

    f = h5.File(file, "r")
    steps_taken = f['steps_taken'][0]
    chain = f['chain'][...][:,:steps_taken:thin,:]
    lnprobability = f['lnprobability'][...][:,:steps_taken:thin]
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

    data = (T, NU, FNU, FERR, UL, INST)
    
    if jetType == 3:
        bounds = fit.boundsCocoon
        fluxFunc = grb.fluxDensityCocoon
    else:
        bounds = fit.boundsJet
        fluxFunc = grb.fluxDensity

    nwalkers = chain.shape[0]
    nsteps = chain.shape[1]
    ndim = chain.shape[2]
    nburn = nsteps/4

    chainBurned = chain[:,nburn:,:]
    lnprobabilityBurned = lnprobability[:,nburn:]

    flatchainBurned = chainBurned.reshape((-1,ndim))
    flatlnprobabilityBurned = lnprobabilityBurned.reshape((-1,))

    w1 = None

    if jetType != 3:
        derivedVals = [{'label': r"$E_{tot}$", 'f':f_Etot},
                        {'label': r"$\theta_v/\theta_c$", 'f':f_thV_o_thC},
                        {'label': r"$\theta_w/\theta_c$", 'f':f_thW_o_thC}]
        wP = weightPlanck(flatchainBurned)
        wS = weightSHoES(flatchainBurned)
        wL = weightLIGO(flatchainBurned)
        analyze(flatchainBurned, flatlnprobabilityBurned, w1, jetType, X0,
                fitVars, labels, data, "em", optimizeMAP,
                derivedVals)
        analyze(flatchainBurned, flatlnprobabilityBurned, wP, jetType, X0,
                fitVars, labels, data, "em+LIGO+Planck", optimizeMAP,
                derivedVals)
        analyze(flatchainBurned, flatlnprobabilityBurned, wS, jetType, X0,
                fitVars, labels, data, "em+LIGO+SHoES", optimizeMAP,
                derivedVals)
        analyze(flatchainBurned, flatlnprobabilityBurned, wL, jetType, X0,
                fitVars, labels, data, "em+LIGO", optimizeMAP,
                derivedVals)
    else:
        derivedVals = [{'label': r"$E_{tot}$", 'f':f_Etot}]
        wE = weightEtot(flatchainBurned, jetType, X0, fitVars)

        analyze(flatchainBurned, flatlnprobabilityBurned, w1, jetType, X0,
                fitVars, labels, data, "em", optimizeMAP,
                derivedVals)
        analyze(flatchainBurned, flatlnprobabilityBurned, wE, jetType, X0,
                fitVars, labels, data, "em+Ecut", optimizeMAP,
                derivedVals)
        analyze(flatchainBurned[wE==1.0], flatlnprobabilityBurned[wE==1.0],
                None, jetType, X0, fitVars, labels, data, "em+Ecut2", 
                optimizeMAP, derivedVals)


