# IMPORT BLOCK
###############################
###############################
import numpy as np
import astropy.units as u
import astropy.cosmology.units as cu
from astropy.cosmology import Planck18


def soft_to_hard(softflux, gamma, E2hard, E1hard, E2soft, E1soft):
    """
    This is a function that converts from soft X-ray flux to hard X-ray flux
    given the energy ranges and spectral index.

    Args:
        softflux (1D array): An N row array containing soft X-ray flux values
        gamma (float): The spectral index
        E2hard (float): Upper bound of hard energy band
        E1hard (float): Lower bound of hard energy band
        E2soft (float): Upper bound of soft energy band
        E1soft (float): Lower bound of hard energy band
    Returns:
        hardflux (1D array): Outputs converted hard fluxes
        A (1D array): The amplitude of the X-ray flux
    """

    g = 2-gamma #exponent is 2-gamma
    conv = (1.6*10**-9) #conversion factor
    E = E2soft**g - E1soft**g #E2-E1
    A = softflux*g/(E*conv) #the amplitude

    hardflux = conv*A*(E2hard**g - E1hard**g)/g #total equation
    return hardflux, A

def find_mbh(fluxes, redshift, gamma):
    """
    This is a function that finds black hole mass given the flux, redshift, and spectral index.

    Args:
        fluxes (1D array): An N row array containing X-ray flux values
        redshift (1D array): An N row array containing redshifts
        gamma (float): The spectral index
    Returns:
        Mbh_AGN (1D array): Outputs estimated black hole mass for each source
    """
    z_AGN = redshift* cu.redshift #convert to redshift units
    dp = z_AGN.to(u.Mpc, cu.redshift_distance(Planck18, kind='luminosity')) #convert z to distance
    dl = dp*3.08568*10**24 #convert Mpc to cm
    kz = (1+z_AGN)**(gamma-2) #k correction
    Lx_hard = 4 * np.pi * dl * dl * fluxes * kz / (1+z_AGN)**2 #convert flux to luminosity
    Lbol = Lx_hard * 26 #convert to bolometric luminosity
    Ledd = Lbol/1 #assume accreting at Eddington limit
    Mbh_AGN = Ledd/((3*10**4)*3.846*10**33) #Eddington Luminosity to Black hole mass
    Mbh_AGN = Mbh_AGN/(u.Mpc * u.Mpc) #get rid of units
    z_AGN = z_AGN/(cu.redshift)
    z_AGN = z_AGN.astype(np.float64)
    Mbh_AGN = Mbh_AGN.astype(np.float64)
    return Mbh_AGN

#read in eFEDS SDSS
efedsflux2, efedsflux2error, z_AGN = np.loadtxt('efeds_sdss.txt', unpack=True)
efedshardfluxes,efedsamps =soft_to_hard(efedsflux2, 1.8, 10, 2, 2.3, 0.2)
logflux = np.log10(efedshardfluxes)


#find black hole mass
zvals = np.arange(np.min(z_AGN), np.max(z_AGN), 0.02) #zgrid
mbh = find_mbh(np.power(10,logflux), zvals,1.8)
