#/usr/bin/env python

import numpy as np
from astropy.table import Table
import astropy.units as u

import sys
sys.path.append("../../QLFs/")
import Shen20

sys.path.append("../")
from Nqso_Shen_pubtools import Nqso

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)


if len(sys.argv)!=2:
    print("Correct use: python",sys.argv[0]," model")
    sys.exit()
model = sys.argv[1]
if model not in ["A","B"]:
    print("Wrong model: ", model)
    sys.exit()


area = 20000.*u.deg**2

#Create the QLF object.
#model = "B"
#model ="A"
qlf = Shen20.QLF(model=model)

mstar_data = Table.read("mstar_z.vandenberk.{}.dat".format(model), format='ascii')

cato = open("Table10_2.{}.txt".format(model),"w")
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
        #print(Nqso(z[k], z[k+1], mi[i], mi[i+1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=-20, area=area))
        #input()
        Nq = Nqso(z[k], z[k+1], mi[i], mi[i+1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=-20, area=area, cosmo=cosmo)
        cato.write("{0:12.0f}".format(Nq))
        print(Nq)
    cato.write("{0:12.0f}\n".format(Nqso(z[0], z[-1], mi[i], mi[i+1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=-20, area=area, cosmo=cosmo)))

cato.write("{0:7s}".format("Total"))
for k in range(len(z[:-1])):
    cato.write("{0:12.0f}".format(Nqso(z[k], z[k+1], mi[0], mi[-1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=-20, area=area, cosmo=cosmo)))
cato.write("{0:12.0f}\n".format(Nqso(z[0], z[-1], mi[0], mi[-1], 'LSSTi', qlf, mstar_data=mstar_data, Mi_lim=-20, area=area, cosmo=cosmo)))

cato.close()
