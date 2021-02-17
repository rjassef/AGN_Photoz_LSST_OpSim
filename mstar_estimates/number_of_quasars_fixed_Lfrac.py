#!/usr/bin/env python

import numpy as np
from scipy.integrate import dblquad
from astropy.cosmology import Planck13 as cosmo
import astropy.units as u

import Shen20

"""
Function we will integrate. This returns dn/dLfrac * dVc/dz, which is the function we want to integrate to obtain the number of quasars per sterradian.

Note that Lfrac = L/Lstar. So we need to multiply dn/dL by Lstar to obtain dn/dLstar

"""
def dndLfrac_dVcdz(Lfrac,z,qlf):
    Lstar = 10.**(qlf.log_Lstar(z))*qlf.Lstar_units
    dvdc = cosmo.differential_comoving_volume(z)
    return (Lstar*dvdc*qlf.dndL(Lfrac,z)).to(1./u.sr).value


#Create the QLF object.
qlf = Shen20.QLF(model="B")

#Find the number of quasars brighter than 0.001 Lstar between z=0.1 and z=6.0 in the whole sky.
zmin = 0.1
zmax = 6.0
Lfrac_min = 1e-3
Lfrac_max = np.inf
N = 4.*np.pi * dblquad(dndLfrac_dVcdz, zmin, zmax, Lfrac_min, Lfrac_max, args={qlf})[0]
print(int(N))
