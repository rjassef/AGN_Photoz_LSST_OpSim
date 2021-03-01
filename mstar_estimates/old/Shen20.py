import numpy as np
from astropy.constants import L_sun
import astropy.units as u

class QLF(object):

    """
    Class that implements the Quasar Luminosity Function from Shen et al. (2020, arXiv:2001.02696), shortenned to S20 hereafter.

    Each of the functional form parameters, as well as the space density, are implemented as methods. Each method is documented with comments throughout its definition.

    Parameters
    ----------
    model: string
        A or B. Determins the Shen et al. QLF parametrization to use. Default is model B.

    """

    def __init__(self, model="B"):

        #Save the model.
        self.model = model

        #Check which model is being used.
        if model=="A":
            icol = 1
        elif model=="B":
            icol = 4
        else:
            print("Wrong model")
            return

        #Read the QLF model parameters from Shen20.dat, which is just Table 4 of S20.
        #self.a0, self.a1, self.a2, self.b0, self.b1, self.b2, self.c0, self.c1, self.c2, self.d0, self.d1 = np.loadtxt("Shen20.dat", usecols=[icol])
        S20_T4 = open("Shen20.dat")
        for line in S20_T4:
            x = line.split()
            exec("self.{0} = {1}".format(x[0], x[icol]))
        S20_T4.close()

        #Reference redshift parameter (see section 4.2 of S20).
        self.z_ref = 2

        #Units of Lstar and phi_star. Since the methods log_Lstar and log_phi_star return the base 10 logarithm of them, we need to maintain the units in these variables such that we can write
        #
        # Lstar = 10.**(self.log_Lstar(z)) * Lstar_units
        #
        # phi_star = 10.**(self.log_phi_star(z)) * phi_star_units
        #
        # to get those parameters in the correct units.
        self.Lstar_units = L_sun
        self.phi_star_units = u.dex**-1 * u.Mpc**-3

        return

    def log_Lstar(self, z):
        """
        This method returns the base 10 logarithm of bolometric value of L_star in units of solar luminosities (L_sun) at redshift z.

        Based on equation (14) of S20.

        """
        x = (1.+z)/(1.+self.z_ref)
        log_Lstar = 2.*self.c0/(x**self.c1+x**self.c2)
        return log_Lstar


    def gamma1(self, z):
        """
        This method returns the value of the gamma_1(z) parameter at redshift z.
        The Chebyshev polynomials are defined later as the method _T(n,x) for n=0, 1 and 2.

        Based on equations (14 and 16) of S20.

        """
        if self.model=="A":
            x = (1+z)
            return self.a0*self._T(0,x) + self.a1*self._T(1,x) + self.a2*self._T(2,x)
        elif self.model=="B":
            x = (1.+z)/(1.+self.z_ref)
            return self.a0 * x**self.a1
        else:
            return np.nan

    def gamma2(self, z):
        """
        This method returns the value of the gamma_2(z) parameter at redshift z.

        Based on equation (14) of S20.

        """
        x = (1.+z)/(1.+self.z_ref)
        return 2*self.b0/(x**self.b1+x**self.b2)

    def log_phi_star(self, z):
        """
        This method returns the base 10 logarithm of phi_star(z) in units of mag^-1 cMpc^-3.

        Based on equation (14) of S20.

        """
        x = (1+z)
        return self.d0*self._T(0,x) + self.d1*self._T(1,x)


    def phi_bol(self, L, z):
        """
        This method returns the space density of quasars with bolometric luminosity L, in units of mag^-1 cMpc^-3.

        Based on equation (11) of S20.

        """
        Lstar = 10.**(self.log_Lstar(z)) * self.Lstar_units
        x = L/Lstar
        phi_star = 10.**(self.log_phi_star(z)) * self.phi_star_units
        return phi_star/(x**self.gamma1(z)+x**self.gamma2(z))

    def dndL(self, Lfrac, z):
        """
        This method returns the number density of quasars with bolometric luminosity L, in units of Lsun^-1 cMpc^-3.

        Based on equation (12) of S20.

        Parameters
        ----------
        Lfrac: float
            Bolometric luminosity fraction L/L*.

        z: float
            Redshift

        """
        Lstar = 10.**(self.log_Lstar(z)) * self.Lstar_units
        phi_star = 10.**(self.log_phi_star(z)) * self.phi_star_units
        phi_star_prime = phi_star/np.log(10.) * u.dex
        alpha = -(self.gamma1(z)+1)
        beta  = -(self.gamma2(z)+1)
        return (phi_star_prime/Lstar)/(Lfrac**(-alpha)+Lfrac**(-beta))

    def dndLfrac(self, Lfrac, z):
        """
        This method returns the number density of quasars with bolometric luminosity fraction Lfrac=L/L*, in units of Lsun^-1 cMpc^-3.

        Based on equation (12) of S20, modified to make the derivative be with respect to L/L* instead of L.

        This is just a simple rewrite of equation (12) to avoid calculating L* in every iteration if we know the value of L/L* instead of L. This is useful for integrating over L/L* instead of over L.

        Parameters
        ----------
        Lfrac: float
            Bolometric luminosity fraction L/L*.

        z: float
            Redshift

        """
        phi_star = 10.**(self.log_phi_star(z)) * self.phi_star_units
        phi_star_prime = phi_star/np.log(10.) * u.dex
        alpha = -(self.gamma1(z)+1)
        beta  = -(self.gamma2(z)+1)
        return phi_star_prime/(Lfrac**(-alpha)+Lfrac**(-beta))


    def L1450(self, Lbol):
        """
        This method returns the L1450 monochromatic luminosity of a quasar of bolometric luminosity Lbol using equation (5) and the coefficients in Table 1 of S20.

        While the typical use of equation (5) is to determine Lbol given an observable monochromatic luminosity, here we use the conversion to go from Lbol to L1450. A direct application of this function is used in the accompanying script mstar.vandenberk.py, where we want to estimate the observed fluxes of a type 1 quasar with bolometric luminosity equal to L*.

        Note that the monochromatic luminosity is defined as in Table 1 of S20, so the units are the same as those in Lbol. In other words, this method return nu*L_nu, not L_nu.

        """
        #Coefficients from Table 1 for the UV luminosity.
        c1, k1, c2, k2 = 1.862, -0.361, 4.870, -0.0063
        #Implementation of equation (5).
        x = Lbol/(1e10*L_sun)
        bc = c1*x**k1 + c2*x**k2
        return Lbol/bc


    def _T(self, n, x):
        """
        This method returns the nth order Chebyshev polynomial Tn(x) evaluated at the input value x, with n=1, 2 or 3.

        See eqn. (14) fo S20.
        """
        if n==0:
            return 1.0
        elif n==1:
            return x
        elif n==2:
            return 2*x**2-1
        else:
            raise ValueError("Chebyshev polynomial not implemented for order n={}. Returning 0.".format(n))
            return 0
