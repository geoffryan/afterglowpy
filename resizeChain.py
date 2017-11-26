import sys
import h5py as h5
import numpy as np

if __name__ == "__main__":

    filename = sys.argv[1]
    nsteps2 = int(sys.argv[2])
    filename2 = sys.argv[3]

    f = h5.File(filename, "a")

    chain = f['chain'][...]
    lnprob = f['lnprobability'][...]
    steps_taken = f['steps_taken'][0]

    keys = []
    vals = []

    for key in list(f):
        if key not in ['chain', 'lnprobability']:
            keys.append(key)
            vals.append(f[key][...])

    f.close()

 
    nwalkers = chain.shape[0]
    nsteps1 = chain.shape[1]
    ndim = chain.shape[2]

    chain2 = np.empty((nwalkers, nsteps2, ndim))
    lnprob2 = np.empty((nwalkers, nsteps2))

    chain2[:,:steps_taken,:] = chain[:,:steps_taken,:]
    lnprob2[:,:steps_taken] = lnprob[:,:steps_taken]
    chain2[:,steps_taken:,:] = 0
    lnprob2[:,steps_taken:] = 0

    f = h5.File(filename2, "w")
    for i,key in enumerate(keys):
        if key not in ['chain', 'lnprobability']:
            f.create_dataset(key, data=vals[i])
    f.create_dataset("chain", data=chain2)
    f.create_dataset("lnprobability", data=lnprob2)
    f.close()

