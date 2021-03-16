#/usr/bin/env python

import numpy as np
from astropy.table import Table
import astropy.units as u

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

import sys
sys.path.append(root_path+"/Quasar_Counts/")
from Fast_Nqso import Fast_Nqso


if len(sys.argv)!=3 and len(sys.argv)!=4:
    print("Correct use: python",sys.argv[0]," qlf_module qlf_model [SED_model]")
    sys.exit()

qlf_module = sys.argv[1]
qlf_model  = sys.argv[2]
if len(sys.argv)==4:
    SED_model = sys.argv[3]
else:
    SED_model = "Richards06"

area = 20000.*u.deg**2

#Create the qso counter object.
qso_counter = Fast_Nqso('LSSTi', qlf_module, qlf_model, SED_model, area)

cato = open("Fast_Table10_2.{0}.{1}.{2}.txt".format(SED_model, qlf_module, qlf_model),"w")
z  = np.array([0.3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.7])
mi = np.array([15.7, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.3])

for i in range(len(mi)):
    cato.write("{0:5.1f} ".format(mi[i]))
cato.write("\n")
for k in range(len(z)):
    cato.write("{0:5.1f} ".format(z[k]))
cato.write("\n")

for i in range(len(mi[:-1])):
    print(mi[i], mi[i+1])
    cato.write("{0:7d}".format(int(mi[i+1])))
    for k in range(len(z[:-1])):
        print(z[k], z[k+1])
        Nq = qso_counter.Nqso(z[k], z[k+1], mi[i], mi[i+1])[0]
        cato.write("{0:12.0f}".format(Nq))
    cato.write("{0:12.0f}\n".format(qso_counter.Nqso(z[0], z[-1], mi[i], mi[i+1])[0]))

cato.write("{0:7s}".format("Total"))
for k in range(len(z[:-1])):
    cato.write("{0:12.0f}".format(qso_counter.Nqso(z[k], z[k+1], mi[0], mi[-1])[0]))
cato.write("{0:12.0f}\n".format(qso_counter.Nqso(z[0], z[-1], mi[0], mi[-1])[0]))

cato.close()
