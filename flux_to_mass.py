# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.cosmology.units as cu
from astropy.cosmology import Planck18
from scipy import interpolate
import colorcet as cc
import matplotlib.image as mpimg
from astropy.coordinates import SkyCoord
import matplotlib.mlab as mlab
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from sklearn.neighbors import KernelDensity
import seaborn as sns
from matplotlib import colors
import astropy
from astropy.stats import freedman_bin_width
from scipy.stats import gaussian_kde
from scipy.stats import rv_histogram
from scipy import stats

import sys
sys.path.insert(0,'./')
from scipy.optimize import fsolve

def soft_to_hard(softflux, gamma, E2hard, E1hard, E2soft, E1soft):
    g = 2-gamma
    conv = (1.6*10**-9)
    E = E2soft**g - E1soft**g
    A = softflux*g/(E*conv)

    hardflux = conv*A*(E2hard**g - E1hard**g)/g
    return hardflux, A
#read in eFEDS SDSS
efedsflux2, efedsflux2error, z_AGN = np.loadtxt('efeds_sdss.txt', unpack=True)
efedshardfluxes,efedsamps =soft_to_hard(efedsflux2, 1.8, 10, 2, 2.3, 0.2)
logflux = np.log10(efedshardfluxes)


#find black hole mass
def find_mbh(fluxes, redshift, gamma):
    z_AGN = redshift* cu.redshift
    dp = z_AGN.to(u.Mpc, cu.redshift_distance(Planck18, kind='luminosity')) #convert z to distance
    dl = dp*3.08568*10**24 #convert Mpc to cm
    kz = (1+z_AGN)**(gamma-2)
    Lx_hard = 4 * np.pi * dl * dl * fluxes * kz / (1+z_AGN)**2
    Lbol = Lx_hard * 26
    Ledd = Lbol/1
    Mbh_AGN = Ledd/((3*10**4)*3.846*10**33)
    Mbh_AGN = Mbh_AGN/(u.Mpc * u.Mpc)
    z_AGN = z_AGN/(cu.redshift)
    z_AGN = z_AGN.astype(np.float64)
    Mbh_AGN = Mbh_AGN.astype(np.float64)
    return Mbh_AGN

zvals = np.arange(np.min(z_AGN), np.max(z_AGN), 0.02) #zgrid
mbh = find_mbh(np.power(10,logflux), zvals,1.8)
