import numpy as np
from astropy.constants import L_sun
import astropy.units as u

class QLF(object):

    """
    Quasar Luminosity Function from Shen et al. (2020, arXiv:2001.02696)

    """

    def __init__(self, model="B"):

        if model=="A":
            icol = 1
        elif model=="B":
            icol = 4
        else:
            print("Wrong model")
            return
        self.a0, self.a1, self.a2, self.b0, self.b1, self.b2, self.c0, self.c1, self.c2, self.d0, self.d1 = np.loadtxt("Shen20.dat", usecols=[icol])

        self.z_ref = 2

        self.Lstar_units = L_sun
        self.phi_star_units = u.dex**-1 * u.Mpc**-3

        return

    def log_Lstar(self, z):
        x = (1.+z)/(1.+self.z_ref)
        log_Lstar = 2.*self.c0/(x**self.c1+x**self.c2)
        return log_Lstar

    def _T(self, n, x):
        if n==0:
            return 1.0
        elif n==1:
            return x
        elif n==2:
            return 2*x**2-1
        else:
            return 0

    def gamma1(self, z):
        x = (1+z)
        return self.a0*self._T(0,x) + self.a1*self._T(1,x) + self.a2*self._T(2,x)

    def gamma2(self, z):
        x = (1.+z)/(1.+self.z_ref)
        return 2*self.b0/(x**self.b1+x**self.b2)

    def log_phi_star(self, z):
        x = (1+z)
        return self.d0*self._T(0,x) + self.d1*self._T(1,x)


    def phi_bol(self, L, z):
        Lstar = 10.**(self.log_Lstar(z)) * self.Lstar_units
        x = L/Lstar
        phi_star = 10.**(self.phi_star(z)) * self.phi_star_units
        return phi_star/(x**self.gamma1(z)+x**self.gamma2(z))


    def L1450(self, Lbol):

        c1, k1, c2, k2 = 1.862, -0.361, 4.870, -0.0063

        x = Lbol/(1e10*L_sun)
        bc = c1*x**k1 + c2*x**k2
        return Lbol/bc
