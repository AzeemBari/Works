# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time
from astropy.io import fits
import sys

sys.path.insert(0, '../../../runnoe/week8/')
import sdss_sweep_data_index
from sdss_sweep_data_index import sdss_sweep_data_index

def positions(filename):
     hdu = fits.open(filename)
     obj = hdu[1].data
     ra = obj["RA"]
     dec = obj["DEC"]
     return ra,dec

def get_qsos(ra,dec,rad):
    swfiles = sdss_sweep_data_index(ra, dec, rad, objtype='star',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final/')
    #AB want primary sources
    for i in range(0, len(swfiles)):
        swfiles[i] = swfiles[i].replace("star","stargal-primary")
    #AB grab the RAs, decs, and fluxes for these objects
    hdus = [fits.open(swfile) for swfile in swfiles]
    objs = [hdu[1].data for hdu in hdus]
    objs = np.hstack(objs)
    sdss_ra = objs["RA"]
    sdss_dec = objs["DEC"]
    flux = objs["PSFFLUX"]
    points = objs["OBJC_TYPE"] #make sure they are point sources!
    flags = objs["FLAGS"]


    #qso_ra, qso_dec = positions("../../../runnoe/week11/qsos-ra180-dec30-rad3.fits")
    qso_ra, qso_dec = positions("../../../runnoe/week13/qsos-ra210-dec48-rad2.fits")
    print(len(qso_ra))
    #AB now for WISE object
    for i in range(0, len(swfiles)):
            swfiles[i] = swfiles[i].replace("star","wise-star")
    wisehdus = [fits.open(swfile) for swfile in swfiles]
    wiseobjs = [hdu[1].data for hdu in wisehdus]
    wiseobjs = np.hstack(wiseobjs)
    #AB grab the fluxes in W1 and W2
    wiseflux = [wiseobjs["W1_NANOMAGGIES"],wiseobjs["W2_NANOMAGGIES"]]



    dtype = np.dtype([('RA', 'f8'), ('DEC', 'f8'), ('PSFFLUX', 'f8', (len(objs["PSFFLUX"][0]))), ("EXTINCTION", 'i8', (len(objs["EXTINCTION"][0]))), ("RESOLVE_STATUS", 'i8'), ("OBJC_TYPE", 'i8'), ("OBJC_FLAGS", 'i8'), ("OBJC_FLAGS2", 'i8'), ("FLAGS", 'i8', (len(objs["FLAGS"][0]))), ("FLAGS2", 'i8', (len(objs["FLAGS2"][0]))), ("W1", 'f8'), ("W2", 'f8')])
    newdat = np.zeros(len(objs), dtype=dtype)

    cols = ["RA", "DEC", "PSFFLUX", "EXTINCTION", "RESOLVE_STATUS", "OBJC_TYPE", "OBJC_FLAGS", "OBJC_FLAGS2", "FLAGS", "FLAGS2"]
    for str in cols:
       newdat[str] = objs[str]

    newdat["W1_NANOMAGGIES"] = wiseobjs["W1_NANOMAGGIES"]
    newdat["W2_NANOMAGGIES"] = wiseobjs["W2_NANOMAGGIES"]

    c1 = SkyCoord(ra=sdss_ra, dec=sdss_dec, frame='icrs', unit='deg')
    c2 = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')
    sep = c2.separation(c1)
    circ = np.where(sep.degree<rad)
    testdata = newdat[circ]
    return(testdata)
if __name__ == '__main__':
    testdata = get_qsos(210,48,2)
