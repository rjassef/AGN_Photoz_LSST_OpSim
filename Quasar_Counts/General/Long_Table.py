#/usr/bin/env python

import numpy as np
from astropy.table import Table
import astropy.units as u
import multiprocessing as mp
from functools import partial

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

import sys
sys.path.append(root_path+"/QLFs/")

sys.path.append(root_path+"/Quasar_Counts/")
from Nqso_v2_MiLim import Nqso

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

sys.path.append(root_path+"/Filter_Curves/")
from lam_eff import lam_eff

###

def get_Nqso(z, m, LSSTfilter, qlf, mstar_data, Mi_lim, area, cosmo, m_index_use):
    N = np.zeros((len(m_index_use),len(z)-1))
    for j,i in enumerate(m_index_use):
        for k in range(len(z[:-1])):
            N[j,k] = Nqso(z[k], z[k+1], m[i], m[i+1], LSSTfilter, qlf, area=area, mstar_data=mstar_data, Mi_lim=Mi_lim, cosmo=cosmo)
    return N

if len(sys.argv)!=4 and len(sys.argv)!=5:
    print("Correct use: python",sys.argv[0]," qlf_module qlf_model filter [SED_model]")
    sys.exit()

qlf_module = sys.argv[1]
qlf_model  = sys.argv[2]
LSSTfilter = sys.argv[3]
if len(sys.argv)==5:
    SED_model = sys.argv[4]
else:
    SED_model = "Richards06"

#Create the QLF object.
exec("import {}".format(qlf_module))
exec("qlf = {0}.QLF(model=\"{1}\")".format(qlf_module, qlf_model))

#Define the area so that we get the number per square degrees.
area = 1.0*u.deg**2

#Read the appropriate mstar data file.
mstar_data = Table.read(root_path+"/mstar_estimates/mstar_z.{0}.{1}.{2}.dat".format(SED_model, qlf_module, qlf_model), format='ascii')

#Only allow the redshifts ranges to get up to the point where the effective wavelength of the filter is longwards of the Lyman break.
dz = 0.1
zmin = 0.1
zmax = np.min([7.0, (lam_eff[LSSTfilter]/(912.*u.AA)).to(1.).value])
z  = np.arange(zmin, zmax+0.1*dz, dz)
#z  = np.arange(0.1, 7.0, 0.1)
m  = np.arange(15.7, 26.3, 0.1)

#z  = np.arange(5.0, 6.0, 0.2)
#z  = np.arange(1.0, 2.0, 0.2)
#m  = np.arange(20.0, 22.0, 0.2)

Ncpu = mp.cpu_count()-1
m_index_use = np.arange(len(m)-1)
m_index_use_split = np.array_split(m_index_use, Ncpu)

Mi_lim = -20
Pool = mp.Pool(Ncpu)
func = partial(get_Nqso, z, m, LSSTfilter, qlf, mstar_data, Mi_lim, area, cosmo)

Output = Pool.map(func, m_index_use_split)
Output = np.vstack(Output)
Pool.close()

cato = open("Long_Table.{0}.{1}.{2}.{3}.txt".format(LSSTfilter, qlf_module, qlf_model, SED_model),"w")
for i in range(len(m)):
    cato.write("{0:5.1f} ".format(m[i]))
cato.write("\n")
for k in range(len(z)):
    cato.write("{0:5.1f} ".format(z[k]))
cato.write("\n")

for i in range(len(m[:-1])):
    for k in range(len(z[:-1])):
        cato.write("{0:12.2f}".format(Output[i,k]))
    cato.write("\n")

cato.close()
