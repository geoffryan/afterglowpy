import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from grbpy import c, mp, me, ee, sigmaT

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
            'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
            'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
            'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

def calcSpec(nu, g, te, n0, p, epse, epsB, ksiN, IC=True):

    u2 = (g-1)*(g+1)
    beta = math.sqrt(u2)/g
    d = math.sqrt((1-beta)/(1+beta))
    a = 1-beta

    n = 4*g*n0
    eth = (g-1) * n * mp*c*c
    B = math.sqrt(8*np.pi*eth*epsB)
    em = ksiN * n * B

    gm = (2.0-p)/(1.0-p) * epse * eth / (ksiN * n * me*c*c)
    gcs = 6*np.pi * g * me*c / (sigmaT * B*B * te)

    if IC:
        gr = gcs/gm
        y = beta*epse/epsB
        if gr <= 1.0 or gr*gr-gr-a <= 0.0:
            #Fast Cooling (with I.C.)
            X = 0.5*(1 + math.sqrt(1+4*y))
        else:
            #Slow Cooling (with I.C.)
            b = y*math.pow(gr, 2-p)
            Xa = 1+b
            Xb = math.pow(b, 1.0/(4.0-p)) + 1.0/(4.0-p)
            s = b*b/(b*b+1.0)
            X0 = Xa * math.pow(Xb/Xa, s)
            #cc = math.log10(b)/(1.0 - 0.68*p + 0.178*p*p)
            #X0 += X0*(-0.078-0.033*p+0.02*p*p) * math.exp(-0.5*cc*cc)
            X = X0
            imax = 10
            for i in range(imax):
                po = math.pow(X, p-2)
                f = X*X - X - b*po
                df = 2*X - 1 - (p-2)*b*po/X
                dX = -f/df
                X += dX
                if math.fabs(dX/X) < 1e-3:
                    break
            print("IC solve took {0:d} iterations".format(i+1))
    else:
        X = 1.0

    print(X)

    gc = gcs/X

    num = 3*gm*gm * ee * B /(4*np.pi * me*c)
    nuc = 3*gc*gc * ee * B /(4*np.pi * me*c)
 
    nup = d * nu
    Fnu = np.empty(nu.shape)

    if num < nuc:
        segA = (nup < num)
        segB = (nup >= num) * (nup < nuc)
        segC = (nup >= nuc)
        Fnu[segA] = np.power(nup[segA]/num, 1.0/3.0)
        Fnu[segB] = np.power(nup[segB]/num, 0.5*(1.0-p))
        Fnu[segC] = math.pow(nuc/num, 0.5*(1-p)) * np.power(nup[segC]/nuc, -0.5*p)
    else:
        segA = (nup < nuc)
        segB = (nup >= nuc) * (nup < num)
        segC = (nup >= num)
        Fnu[segA] = np.power(nup[segA]/nuc, 1.0/3.0)
        Fnu[segB] = np.power(nup[segB]/nuc, -0.5)
        Fnu[segC] = math.pow(num/nuc, -0.5) * np.power(nup[segC]/num, -0.5*p)

    Fnu *= 2*np.pi * em/(g*g*a*a)

    return Fnu

def testICSolve():

    a = np.logspace(-4, 8, num=101, base=10.0)
    p = np.linspace(2.05, 2.95, 5)
    gcsogm = np.logspace(-6, 6, num=11, base=10.0)

    xp1F = 0.5*(1 + np.sqrt(1+4*a))

    b = a[None,None,:] * np.power(gcsogm[None,:], (2-p)[:,None])[:,:,None]


    cut = 1.0
    xp1S = np.empty(b.shape)
    #xp1S[:,:,:] = np.maximum(1.0, gcsogm)[None,:,None]
    #xp1S[b<=cut] = 1 + (b/(1-b*(p-2)[:,None,None]))[b<=cut]
    #xp1S[b>cut] = np.power(b, 1.0/(4.0-p)[:,None,None])[b>cut]
    #X1 = 1 + (b/(1-b*(p-2)[:,None,None]))

    X1 = 1 + b
    X2 = np.power(b, 1.0/(4.0-p)[:,None,None]) + 1.0/(4.0-p)[:,None,None]

    def switch(x, n):
        return np.power(x,n)/(np.power(x,n) + math.pow(cut,n))
    f = switch(b, 2)
    X0 = X1 * np.power(X2/X1, f)
    cc = np.log10(b)/(1.0 - 0.68*p + 0.178*p*p)[:,None,None]
    X0 += X0*(-0.078-0.033*p+0.02*p*p)[:,None,None] * np.exp(-0.5*cc*cc)
    #X0 = np.power(X1, 1-f) * np.power(X2,f)

    xp1S[:,:,:] = X0[:,:,:]

    imax = 10
    for i in range(imax):
        po = np.power(xp1S, (p-2)[:,None,None])
        f = xp1S*xp1S - xp1S - b*po
        df = 2*xp1S - 1 - b*(p-2)[:,None,None]*po/xp1S
        print(i, f.min(), f.max())
        indmin = np.unravel_index(f.argmin(), b.shape)
        indmax = np.unravel_index(f.argmax(), b.shape)
        print(xp1S[indmin], b[indmin], p[indmin[0]], gcsogm[indmin[1]])
        print(xp1S[indmax], b[indmax], p[indmax[0]], gcsogm[indmax[1]])

        dxox = (f/df) / xp1S
        print(dxox.min(), dxox.max())
        indmin = np.unravel_index(dxox.argmin(), b.shape)
        indmax = np.unravel_index(dxox.argmax(), b.shape)
        print(xp1S[indmin], b[indmin], p[indmin[0]], gcsogm[indmin[1]])
        print(xp1S[indmax], b[indmax], p[indmax[0]], gcsogm[indmax[1]])
        print("")

        xp1S -= f/df

    f = xp1S*xp1S - xp1S - b*np.power(xp1S, (p-2)[:,None,None])
    print(imax, f.min(), f.max())
    f = xp1S*(xp1S - 1 - b*np.power(xp1S, (p-3)[:,None,None]))
    print(imax, f.min(), f.max())
    f = xp1S*(xp1S - 1) - b*np.power(xp1S, (p-2)[:,None,None])
    print(imax, f.min(), f.max())

    imax = np.unravel_index(np.fabs((xp1S-X0)/xp1S).argmax(), X0.shape)
    print(xp1S[imax], X0[imax], (xp1S[imax]-X0[imax])/xp1S[imax])

    X = np.logspace(-6,6,100)
    for i,gr in enumerate(gcsogm):
        fig, ax = plt.subplots(1,1)

        for j,pp in enumerate(p):
            ax.plot(a, xp1S[j,i,:], color=colors[j], lw=2)
            ax.plot(a, np.power(b[j,i,:], 1.0/(4.0-p[j])), ls='-.', color=colors[j], lw=1)
            ax.plot(a, 1+b[j,i,:], ls=':', color=colors[j], lw=1)
            ax.plot(a, X0[j,i,:], ls='-', lw=1, color='k')
            ax.axvline(cut * math.pow(gr, pp-2), lw=1, color=colors[j], ls='--')
        ax.plot(a, xp1F, color='k', lw=2)

        ax.axhline(gr, ls='--', color='k', lw=2)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$\beta \epsilon_e/\epsilon_B$')
        ax.set_ylabel(r'$1+x$')
        fig.savefig("ic_func_{0:d}.png".format(i), dpi=200)
        plt.close(fig)

    err = (xp1S-X0)/(X0)
    errL = xp1S-X0

    def f_errA(pars):
        A = pars[0] + p*pars[1] + p*p*pars[2]
        sig = pars[3] + p*pars[4] + p*p*pars[5]
        cc = np.log10(b)/sig[:,None,None]
        return A[:,None,None] * np.exp(-0.5*cc*cc) * X0

    def chi2(pars):
        myerr = f_errA(pars)
        dif = myerr-errL
        return (dif*dif).sum()

    par0 = np.array([-0.078, -0.033, 0.020, 1.0, -0.68, 0.178])
    #res = opt.minimize(chi2, par0)
    #par = res.x
    #print(res)
    errA = f_errA(par0)

    eAa = np.empty(err.shape)
    eAa = np.power(b, 1.86)
    eAb = np.empty(err.shape)
    eAb[:,:,:] = 1.0/(4.0-p)[:,None,None]
    #errA = np.zeros(b.shape)
    #errA = eAa*eAb / (eAa + eAb)
    #errA = (-0.065 + 0.06*(p-2))[:,None,None] * np.exp(-0.5*(np.log10(b)/0.5)**2) * X0

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, errL[i,:,:].flat, marker='+', ls='', color=colors[i])
        ax.plot(b[i,:,:].flat, errA[i,:,:].flat, marker=',', ls='', color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_d_log.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, errL[i,:,:].flat, marker='+', ls='', color=colors[i])
        ax.plot(b[i,:,:].flat, errA[i,:,:].flat, marker=',', ls='', color='k')
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_d_lin.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs(errL-errA)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_d_corr_log.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, (errL-errA)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_d_corr_lin.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, err[i,:,:].flat, marker='+', ls='', color=colors[i])
        ax.plot(b[i,:,:].flat, (errA/X0)[i,:,:].flat, ls='', marker=',', color='k')
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_lin.png")
    plt.close(fig)
    
    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, err[i,:,:].flat, marker='+', ls='', color=colors[i])
        ax.plot(b[i,:,:].flat, (errA/X0)[i,:,:].flat, ls='', marker=',', color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    #ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err_log.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, err[i,:,:].flat, marker='+', ls='', color=colors[i])
        ax.plot(b[i,:,:].flat, (errA/X0)[i,:,:].flat, ls='', marker=',', color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    ax.set_xlim(1.0e-2, 1.0e2)
    ax.set_ylim(1.0e-4, 1.0e0)
    fig.savefig("ic_err_log_zoom.png")
    plt.close(fig)
    
    derr = np.log(errL[:,:,1:]/errL[:,:,:-1]) / np.log(b[:,:,1:]/b[:,:,:-1])
    bc = np.sqrt(b[:,:,1:]*b[:,:,:-1])

    s1 = np.empty(p.shape)
    s2 = np.empty(p.shape)

    for i,pp in enumerate(p):
        s1[i] = derr[i][(bc[i]>1.0e-5)*(bc[i]<1.0e-2)].mean()
        s2[i] = errL[i][(b[i]>1.0e5)*(b[i]<1.0e8)].mean()
        print(pp, s1[i], s2[i])

    res = np.polyfit(p, s2, 2)
    print(res)
    fig, ax = plt.subplots(1,1)
    ax.plot(p, s2, ls='', marker='+')
    ax.plot(p, p*p*res[0] + p*res[1] + res[2])
    ax.plot(p, -1.0/(p-4))
    fig.savefig("ic_err_s2.png")
    #print(derr[(b>1e-7)*(b<1e-5)].mean(axis=(1,2)))
    #print(derr[(b>1e5)*(b<1e8)].mean(axis=(1,2)))

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(bc[i,:,:].flat, derr[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('linear')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/X0')
    fig.savefig("ic_err_slope.png")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X0-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err0.png")
    plt.close(fig)
        
    po = np.power(X0, (p-2)[:,None,None])
    f = X0*X0 - X0 - b*po
    df = 2*X0 - 1 - b*(p-2)[:,None,None]*po/X0
    X1 = X0 - f/df

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X1-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X0)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err1.png")
    plt.close(fig)
        
    po = np.power(X1, (p-2)[:,None,None])
    f = X1*X1 - X1 - b*po
    df = 2*X1 - 1 - b*(p-2)[:,None,None]*po/X1
    X2 = X1 - f/df

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X2-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X1)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err2.png")
    plt.close(fig)
        
    po = np.power(X2, (p-2)[:,None,None])
    f = X2*X2 - X2 - b*po
    df = 2*X2 - 1 - b*(p-2)[:,None,None]*po/X2
    X3 = X2 - f/df

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X3-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X2)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err3.png")
    plt.close(fig)
        
    po = np.power(X3, (p-2)[:,None,None])
    f = X3*X3 - X3 - b*po
    df = 2*X3 - 1 - b*(p-2)[:,None,None]*po/X3
    X4 = X3 - f/df

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X4-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X2)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err4.png")
    plt.close(fig)
        
    po = np.power(X4, (p-2)[:,None,None])
    f = X4*X4 - X4 - b*po
    df = 2*X4 - 1 - b*(p-2)[:,None,None]*po/X4
    X5 = X4 - f/df

    fig, ax = plt.subplots(1,1)
    for i, pp in enumerate(p):
        ax.plot(b[i,:,:].flat, np.fabs((X5-xp1S)/xp1S)[i,:,:].flat, marker='+', ls='', color=colors[i])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$b$')
    ax.set_ylabel('(x1pS-X2)/x1pS')
    ax.set_ylim(1.0e-16, 1.0e0)
    fig.savefig("ic_err5.png")
    plt.close(fig)
        

if __name__ == "__main__":


    testICSolve()

    nu = np.logspace(6, 20, 200)

    g = 10
    te = 1.0e9
    n0 = 1.0e-3
    p = 2.15
    epse = 0.1
    epsB = 0.001
    ksiN = 1.0


    tes = [1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10]
    ns = [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
    epses = [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
    epsBs = [1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]

    fig, ax = plt.subplots(1,1)

    FnuY = calcSpec(nu, g, te, n0, p, epse, epsB, ksiN, True)
    FnuN = calcSpec(nu, g, te, n0, p, epse, epsB, ksiN, False)
    ax.plot(nu, FnuN, label='Inverse Compton OFF', color=colors[0])
    ax.plot(nu, FnuY, label='Inverse Compton ON', color=colors[1])

    plt.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\nu$")
    ax.set_ylabel(r"$F_\nu$")

    fig.savefig("ic_spec_comp.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(1,1)

    for i,te in enumerate(tes):
    #for i,n0 in enumerate(ns):
    #for i,epse in enumerate(epses):
    #for i,epsB in enumerate(epsBs):
        FnuY = calcSpec(nu, g, te, n0, p, epse, epsB, ksiN, True)
        FnuN = calcSpec(nu, g, te, n0, p, epse, epsB, ksiN, False)
        ax.plot(nu, FnuN, label='Inverse Compton OFF', color=colors[i], ls='--')
        ax.plot(nu, FnuY, label='Inverse Compton ON', color=colors[i], ls='-')

    #plt.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\nu$")
    ax.set_ylabel(r"$F_\nu$")

    fig.savefig("ic_spec.pdf")

