#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.cosmology import Planck13 as cosmo
from astropy.table import Table

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

SED_path     = root_path+"/SEDs/"
Filters_path = root_path+"/Filter_Curves/"

import lrt
import SED_Model

import sys
sys.path.append(root_path+"/QLFs/")
###

"""
This is a bit of a hack that adds the function Lnu to the the lrt_model class. I did it to not modify the published code at this point.

"""

class lrt_model2(SED_Model.lrt_model):

    def Lnu(self, lam):
        if self.comp is None:
            print("Must fit or provide comp before calculating ahat")
            return
        #Check some redshift is provided.
        self.z = self.zspec
        if self.z==None or self.z<0.:
            self.z = self.zspec
            if self.zphot==None:
                print("Must provide a redshift or run pz_fit.")
                return
            else:
                self.z = self.zphot
        #Check kc has been initialized.
        if self._kcinit == False:
            lrt.kcinit("bandmag.dat",1,1,1,self.iverbose)
            self.nchan = lrt.data1b.nchan
            self._kcinit = True
            self._pzinit = False
        #Calculate L5100. Note this is lam*L_lam
        return lrt.get_lnu(self.comp,self.z,lam)


"""
This script creates the table mstar_z.lrt.qlf_module.qlf_model.dat, which holds the apparent magnitude of an L* quasar as a function of redshift for all the LSST bands, as well as the absolute magnitude at 1450A, M_1450, and at i-band.

All output magnitudes are in the AB system.

Command line arguments:
-----------------------
qlf_module: str
    Can be Hopkins07, Shen20, or any mother module in the QLFs folder.

model: str
    QLF model to use (mostly Full for Hopkins07, and A or B for Shen20).

"""

if len(sys.argv)!=3:
    print("Correct use: python",sys.argv[0],"qlf_module model")
    sys.exit()
qlf_module = sys.argv[1]
qlf_model  = sys.argv[2]

exec("import {}".format(qlf_module))

#Create the QLF object.
exec("qlf = {0}.QLF(model=\"{1}\")".format(qlf_module, qlf_model))

#Setup the cosmology file.
cosmo_file = open("cosmo.dat","w")
cosmo_file.write("{}\n".format(cosmo.Om0)) #Omega_Matter
cosmo_file.write("{}\n".format(1.0-cosmo.Om0)) #Omega_Lambda
cosmo_file.write("{}\n".format(0.0)) #Omega_k
cosmo_file.write("{}\n".format(cosmo.H0.value)) #H0
cosmo_file.close()

if qlf_module == "Shen20":
    lam_cal = 1450.*u.AA
    qlf_Lcal = qlf.L1450

elif qlf_module == "Hopkins07":
    lam_cal = 4380.*u.AA
    qlf_Lcal = qlf.L_B

#Read the photometry zero point.
jyzero = np.loadtxt("bandmag.dat",usecols=[2])
filters = np.genfromtxt("bandmag.dat",usecols=[0],dtype='U')

#Get m* for a range in z and all LSST bands.
zmin = 0.01
zmax = 6.0
dlogz = 0.01
zs = 10.**(np.arange(np.log10(1+zmin), np.log10(1+zmax), dlogz)) - 1
mstar = np.zeros((len(zs),len(filters)))
Lnu1450_unnorm = None
for k,z in enumerate(zs):

    #Create the qso object.
    qso = lrt_model2()
    qso.zspec = z
    qso.comp = np.zeros(4)
    qso.comp[0] = 1.0
    qso.ebv = 0.0
    qso.igm = 1.0

    #Get the qso Lnu_1450
    #Lnu1450_unnorm = qso.Lnu(0.145)*u.erg/u.s/u.Hz
    Lnu_cal_unnorm = qso.Lnu(lam_cal.to(u.micron).value)*u.erg/u.s/u.Hz

    #Get L1450 for an L* quasar.
    Lbol = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units
    nu_cal = c/(lam_cal)
    Lnu_cal = (qlf_Lcal(Lbol)/nu_cal).to(u.erg/u.s/u.Hz)

    #Normalize the component luminosities.
    qso.comp *= (Lnu_cal/Lnu_cal_unnorm).to(1.).value

    #Get the unnormalized observed magnitudes.
    qso.get_model_fluxes()
    mstar[k] = -2.5*np.log10(qso.jymod/jyzero)

    #Do not provide magnitudes for bands that potentially do not overlap with the templates.
    cond = lrt.cal1.lbar[:len(filters)]/(1.+qso.zspec) < 0.09
    mstar[k, cond] = 99.


output = np.zeros((len(zs),len(filters)+1))
output[:,0] = zs
output[:,1:]  = mstar
np.savetxt("mstar_z.lrt.{0}.{1}.dat".format(qlf_module, qlf_model), output,fmt='%15.5f')

import matplotlib.pyplot as plt
for j,filter in enumerate(filters):
    cond = mstar[:,j]<99.
    plt.plot(zs[cond],mstar[cond,j],label=filter)
plt.legend()
plt.xscale('log')
plt.ylim([13.,25.])
plt.xlabel('Redshift')
plt.ylabel('Observed magnitude of L* quasar (AB)')
plt.title('LRT AGN template, {0} QLF, {1} model'.format(qlf_module, qlf_model))
plt.savefig("mstar_z.lrt.{0}.{1}.png".format(qlf_module, qlf_model))
