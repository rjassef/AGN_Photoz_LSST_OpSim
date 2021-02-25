import numpy as np
from scipy.interpolate import interp2d
import astropy.units as u

####

#All this gets executed at the time the code is imported.

cat = open("Long_Table.txt")
mi = np.array([float(x) for x in cat.readline().split()])
zs = np.array([float(x) for x in cat.readline().split()])
miz_data = np.loadtxt(cat)
cat.close()

#Make the table cumulative.
c_miz_data = np.zeros((miz_data.shape[0]+1, miz_data.shape[1]+1))
c_miz_data[1:,1:] = miz_data
c_miz_data = np.cumsum(c_miz_data, axis=0)
c_miz_data = np.cumsum(c_miz_data, axis=1)

#Now, make a 2D interpolation.
Nqso_cumulative = interp2d(zs, mi, c_miz_data, kind='cubic')

#Define the function that returns the number of quasars.
def Nqso(zmin, zmax, mi_bright, mi_faint, filter, area=4.*np.pi*u.sr):

    if filter!='LSSTi':
        print("Only i-band is working for the moment.")
        return np.nan

    N11 = Nqso_cumulative(zmin, mi_bright)[0]
    N12 = Nqso_cumulative(zmax, mi_bright)[0]
    N21 = Nqso_cumulative(zmin, mi_faint )[0]
    N22 = Nqso_cumulative(zmax, mi_faint )[0]

    return (N22 - N21 - N12 + N11)*(area/(4.*np.pi*u.sr)).to(1)
