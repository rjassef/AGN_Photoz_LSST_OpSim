#!/usr/bin/env python

from astropy.table import Table
import astropy.units as u
from scipy.interpolate import interp1d

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

class P92_Extinction(object):

    def __init__(self, red_type):

        self.red_type = red_type
        if red_type not in ["MW", "LMC", "SMC"]:
            print("Model {} not recognized.".format(red_type))
            return

        #Start by reading the RV tables.
        self.R_V = dict()
        RV_file = open(root_path+"/QLFs/Pei92_RV.txt")
        for line in RV_file:
            x = line.split()
            self.R_V[x[0]] = float(x[1])
        RV_file.close()

        #Now, read the extinction tables.
        extinction_curves = Table.read(root_path+"/QLFs/Pei92_Extinction_curves.txt", format='ascii')

        #For the requested model, caculate Xi and make the interpolation function.
        lam     = 1./extinction_curves['ilam_'+red_type] * u.um
        Erat    = extinction_curves['E_rat_'+red_type]
        self.xi_func = interp1d(lam.value, (Erat+self.R_V[red_type])/(1+self.R_V[red_type]), fill_value='extrapolate')

        return

    def xi(self,lam):
        return self.xi_func(lam.to(u.um).value)
