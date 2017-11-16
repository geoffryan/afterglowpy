import numpy as np
import json
import requests

h = 6.62607004e-27
eV = 1.60218e-12

keV = 1.0e3 * eV

class OKCData:

    alldata = None
    name = ""
    t0 = 57982.528519 #MJD
    
    def __init__(self, objname):
        self.name = objname
        self.filename = objname + ".json"

        try:
            f = open(self.filename, "r")
            print("Loading from: " + self.filename)
            res_dict = json.load(f)
            f.close()
            self.alldata = res_dict[objname]['photometry']
        except:
            print("Requesting data from Open Kilonova Catalog (OKC)")
            res = requests.get(
            "https://api.kilonova.space/{0:s}/photometry/".format(
                objname))
            res_dict = res.json()
            f = open(self.filename, "w")
            json.dump(res_dict, f)
            f.close()
            self.alldata = res_dict[objname]['photometry']


    def filterData(self, tok):

        result = []

        for obj in self.alldata:
            if tok in obj.keys():
                result.append(obj)

        return result

    def getRadio(self):

        radioList = self.filterData('fluxdensity')

        N = len(radioList)

        t = np.empty(N, dtype=np.float)
        nu = np.empty(N, dtype=np.float)
        Fnu = np.empty(N, dtype=np.float)
        eFnu = np.empty(N, dtype=np.float)

        for i,obj in enumerate(radioList):
            time = float(obj['time'])
            freq = float(obj['frequency'])
            flux = float(obj['fluxdensity'])

            if obj['source'] == '20':
                # Alexander et al.
                if 'upperlimit' in obj.keys():
                    flux = 0.0
                    eflux = (float(obj['e_fluxdensity'])
                                / float(obj['upperlimitsigma']))
                else:
                    eflux = float(obj['e_fluxdensity'])
            else:
                if 'upperlimit' in obj.keys():
                    flux = 0.0
                    eflux = (float(obj['fluxdensity'])
                                / float(obj['upperlimitsigma']))
                else:
                    eflux = float(obj['e_fluxdensity'])

            if obj['u_fluxdensity'] == u'\xb5Jy':
                #Catalog was in micro-Jy
                flux *= 1.0e-3
                eflux *= 1.0e-3

            if obj['u_time'] == 'MJD':
                time = 86400.0 * (time - self.t0)

            if obj['u_frequency'] == 'GHz':
                freq *= 1.0e9

            t[i] = time
            nu[i] = freq
            Fnu[i] = flux
            eFnu[i] = eflux
        
        return t, nu, Fnu, eFnu

    def getXRay(self):

        xrayList = self.filterData('energy')

        N = len(xrayList)

        t = np.empty(N, dtype=np.float)
        nu = np.empty(N, dtype=np.float)
        Fnu = np.empty(N, dtype=np.float)
        eFnu = np.empty(N, dtype=np.float)

        for i,obj in enumerate(xrayList):
            time = np.array([float(x) for x in obj['time']])
            en = np.array([float(x) for x in obj['energy']])
            flux = float(obj['unabsorbedflux'])

            if 'upperlimit' in obj.keys():
                #Assuming 3-sigma upper limits
                flux = 0.0
                eflux = np.array([(float(obj['unabsorbedflux']))])/3.0
            else:
                eflux = np.array([float(obj['e_lower_unabsorbedflux']),
                            float(obj['e_upper_unabsorbedflux'])])

            if obj['u_time'] == 'MJD':
                time = 86400.0 * (time - self.t0)

            if obj['u_energy'] == 'keV':
                en *= keV
            else:
                en *= eV

            if obj['u_flux'] == 'ergs/s/cm^2':
                flux *= 1.0e23
                eflux *= 1.0e23

            freq = en/h

            t[i] = time.mean()
            nu[i] = freq.mean()
            Fnu[i] = flux / (freq[1] - freq[0])
            eFnu[i] = eflux.mean() / (freq[1] - freq[0])
        
        return t, nu, Fnu, eFnu

def getUnique(lis):

    u = [lis[0]]
    for a in lis:
        if a not in u:
            u.append(a)
    return u

if __name__ == "__main__":

    dat = OKCData("GW170817")

    tR, nuR, FnuR, eFnuR = dat.getRadio()
    tX, nuX, FnuX, eFnuX = dat.getXRay()

    print(len(tR))
    print(len(tX))

    import matplotlib.pyplot as plt

    realR = FnuR>0
    limitR = FnuR==0.0
    realX = FnuX>0
    limitX = FnuX==0.0

    fig, ax = plt.subplots(1,2)
    ax[0].errorbar(tR[realR], FnuR[realR], eFnuR[realR], ls='')
    ax[0].plot(tR[limitR], eFnuR[limitR], ls='', marker='v')
    ax[0].errorbar(tX[realX], FnuX[realX], eFnuX[realX], ls='')
    ax[0].plot(tX[limitX], eFnuX[limitX], ls='', marker='v')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[1].errorbar(nuR[realR], FnuR[realR], eFnuR[realR], ls='')
    ax[1].plot(nuR[limitR], eFnuR[limitR], ls='', marker='v')
    ax[1].errorbar(nuX[realX], FnuX[realX], eFnuX[realX], ls='')
    ax[1].plot(nuX[limitX], eFnuX[limitX], ls='', marker='v')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    plt.show()
