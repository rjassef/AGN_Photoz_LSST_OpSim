import numpy as np
from astropy.constants import L_sun
import astropy.units as u

import os
import re
root_path = re.search("(.*/AGN_Photoz_LSST_OpSim)/*",os.getcwd()).group(1)

import sys
sys.path.append(root_path+"/QLFs/")
from Pei92 import P92_Extinction

class QLF(object):

    """
    Class that implements the Quasar Luminosity Function from Hopkins et al. (2007, ApJ, 654. 731), shortenned to H07 hereafter.

    Each of the functional form parameters, as well as the space density, are implemented as methods. Each method is documented with comments throughout its definition.

    Parameters
    ----------
    model: string
        LDDE, PDE or Full. Determins the H07 QLF parametrization to use. Default is Full.

    """

    def __init__(self, model="LDDE"):

        #Save the model call.
        self.model = model

        #Check which model is being used.
        if model=="PDE":
            icol = 1
        elif model=="LDDE":
            icol = 3
        elif model=="Full":
            icol = 5
        else:
            print("Wrong model")
            return

        #Read the QLF model parameters from Shen20.dat, which is just Table 4 of H07.
        H07_T4 = open(root_path+"/QLFs/Hopkins07.dat")
        for line in H07_T4:
            x = line.split()
            exec("self.{0} = {1}".format(x[0],x[icol]))
        H07_T4.close()

        #Reference redshift parameter (see section 4.2 of S20).
        self.z_ref = 2

        #Dust to gas ratio assumed: A_B/NH
        self.dgr_local = 8.47e-22 * u.cm**2

        #Units of Lstar and phi_star. Since the methods log_Lstar and log_phi_star return the base 10 logarithm of them, we need to maintain the units in these variables such that we can write
        #
        # Lstar = 10.**(self.log_Lstar(z)) * Lstar_units
        #
        # phi_star = 10.**(self.log_phi_star(z)) * phi_star_units
        #
        # to get those parameters in the correct units.
        self.Lstar_units = L_sun
        self.phi_star_units = u.dex**-1 * u.Mpc**-3

        #The table has log_Lc, but it is more useful to have Lc.
        self.Lc = 10.**self.log_Lc * self.Lstar_units

        #Set the reddening model.
        self.red_model = P92_Extinction("MW")

        #Coefficients to calculate the bolometric correction for B-band.
        self.c_B = np.array([6.25, 9.00])
        self.k_B = np.array([-0.37, -0.012])

        #Coefficients to calculate the bolometric correction dispersion for B-band.
        self.sig1_B, self.beta_B, self.sig2_B = 0.08, -0.25, 0.060

        return

    def log_Lstar(self, z):
        """
        This method returns the base 10 logarithm of bolometric value of L_star in units of solar luminosities (L_sun) at redshift z.

        Based on equation (9) of H07.

        """
        x = (1.+z)/(1.+self.z_ref)
        xi = np.log10(x)
        log_Lstar = self.log_Lstar0 + self.k_L1*xi + self.k_L2*xi**2 + self.k_L3*xi**3
        return log_Lstar


    def gamma1(self, z):
        """
        This method returns the value of the gamma_1(z) parameter at redshift z.

        Based on equation (17) of H07.

        """
        x = (1+z)/(1.+self.z_ref)
        return self.gamma1_0 * x**(self.k_gamma1)

    def gamma2(self, z):
        """
        This method returns the value of the gamma_2(z) parameter at redshift z.

        Based on equation (20) of H07.

        """
        x = (1.+z)/(1.+self.z_ref)
        return 2*self.gamma2_0/(x**self.k_gamma2_1+x**self.k_gamma2_2)


    def log_phi_star(self, z):
        """
        This method returns the base 10 logarithm of phi_star(z) in units of mag^-1 cMpc^-3.

        In H07 this is a constant, so this is just for compatibility whith the routines that also use the Shen20 QLF.

        """
        if isinstance(z, np.ndarray):
            return self.log_phistar*np.ones(len(z))
        else:
            return self.log_phistar


    def e_d(self, L, z):
        """
        This method returns the value of e_d(L,z). Only useful for the PDE and LDDE models.

        Based on equations (13-16) of H07.

        """
        z_c = z_c0
        if L<=self.Lc:
            z_c *= ((L/self.Lc).to(1))**self.alpha

        p1 = self.p1_46 + self.beta_1*np.log10(L/(1e46*u.erg/u.s))
        p2 = self.p2_46 + self.beta_2*np.log10(L/(1e46*u.erg/u.s))

        if z<z_c:
            return (1.+z)**self.p1
        else:
            return (1.+z)**self.p1 * ((1.+z)/(1.+self.z_ref))**p2


    def phi_bol(self, L, z):
        """
        This method returns the space density of quasars with bolometric luminosity L, in units of dex^-1 cMpc^-3.

        Based on equation (11) of S20.

        """
        Lstar = 10.**(self.log_Lstar(z)) * self.Lstar_units
        x = L/Lstar
        phi_star = 10.**(self.log_phi_star(z)) * self.phi_star_units
        return phi_star/(x**self.gamma1(z)+x**self.gamma2(z))

    def phi_bol_Lfrac(self, Lfrac, z):
        """
        This method returns the space density of quasars with bolometric luminosity L, in units of dex^-1 cMpc^-3.

        Based on equation (11) of S20.

        """
        phi_star = 10.**(self.log_phi_star(z)) * self.phi_star_units
        return phi_star/(Lfrac**self.gamma1(z)+Lfrac**self.gamma2(z))


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


    def L_B(self, Lbol):
        """
        This method returns the B-band luminosity of a quasar of bolometric luminosity Lbol using equation (2) of H07.

        While the typical use of equation (2) is to determine Lbol given an observable monochromatic luminosity, here we use the conversion to go from Lbol to LB. A direct application of this function is used in the accompanying script mstar.vandenberk.py, where we want to estimate the observed fluxes of a type 1 quasar with bolometric luminosity equal to L*.

        """
        #Coefficients from Table 1 for the UV luminosity.
        c1, k1, c2, k2 = 6.25, -0.37, 9.00, -0.012
        #Implementation of equation (5).
        x = Lbol/(1e10*L_sun)
        bc = c1*x**k1 + c2*x**k2
        return Lbol/bc

    def L_x(self, Lbol):
        """
        This method returns the hard X-rayluminosity of a quasar of bolometric luminosity Lbol using equation (2) of H07.

        """
        #Coefficients from Table 1 for the UV luminosity.
        c1, k1, c2, k2 = 10.83, 0.28, 6.08, -0.020
        #Implementation of equation (5).
        x = Lbol/(1e10*L_sun)
        bc = c1*x**k1 + c2*x**k2
        return Lbol/bc

    def L_x_Lfrac(self, Lfrac, Lstar_10):
        """
        This method returns the hard X-rayluminosity of a quasar of bolometric luminosity Lbol using equation (2) of H07.

        """
        #Coefficients from Table 1 for the UV luminosity.
        c1, k1, c2, k2 = 10.83, 0.28, 6.08, -0.020
        #Implementation of equation (5).
        x = Lfrac * Lstar_10
        bc = c1*x**k1 + c2*x**k2
        return Lfrac*Lstar_10*1e10*L_sun/bc


    def xi(self, lam):
        return self.red_model.xi(lam)


    def fNH(self, log_NH, Lfrac, Lstar_10, z):

        #Get the hard x-ray luminosity for each Lfrac in units of 10^44 erg/s. This will be useful later.
        Lx = self.L_x_Lfrac(Lfrac, Lstar_10)
        lLx_44 = np.log10(Lx/(1e44*u.erg/u.s)).value

        #logLx = np.log10(self.L_x(Lfrac*10.**(self.log_Lstar(z))*self.Lstar_units)/(1e44*u.erg/u.s))
        eps = 1.7
        psi_max = (1.+eps)/(3+eps)
        psi_44  = 0.47
        beta_L  = 0.1
        psi = psi_44 - beta_L*lLx_44
        psi = np.where(psi<0., 0., psi)
        psi = np.where(psi>psi_max, psi_max, psi)

        f_low = 2.0 - ((5.+2.*eps)/(1.+eps))*psi
        f_med = (1./(1.+eps))*psi
        f_hig = (eps/(1.+eps))*psi
        f_compton = f_hig

        f_low = f_low / (1. + f_compton)
        f_med = f_med / (1. + f_compton)
        f_hig = f_hig / (1. + f_compton)

        f_NH = np.zeros((len(log_NH), len(lLx_44)))
        f_low = np.tile(f_low, [len(log_NH),1])
        f_med = np.tile(f_med, [len(log_NH),1])
        f_hig = np.tile(f_hig, [len(log_NH),1])
        f_compton = np.tile(f_compton, [len(log_NH),1])

        log_NH_2D = np.tile(log_NH, [len(lLx_44),1]).T

        f_NH = np.where((log_NH_2D>=20.0) & (log_NH_2D<=20.5), f_low, f_NH)
        f_NH = np.where((log_NH_2D>20.5) & (log_NH_2D<=23.0), f_med, f_NH)
        f_NH = np.where((log_NH_2D>23.0) & (log_NH_2D<=24.0), f_hig, f_NH)
        f_NH = np.where((log_NH_2D>24.0) & (log_NH_2D<=25.0), f_compton, f_NH)

        return f_NH

    def get_sigma(self, Lfrac, Lstar_10):
        """
        This function calculates the dispersion of the bolometric correction for B-band.

        Parameters
        ----------

        Lfrac: numpy array
            Values of L/Lstar for which to calculate Lfrac_lam = L_lam/L_lam(Lstar)

        Lstar_10: float
            Value of Lstar in units of 10^10 Lsun.

        """
        return self.sig1_B*(Lstar_10*10 * Lfrac)**self.beta_B + self.sig2_B

    def dgr(self, z):
        return self.dgr_local
