# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import time


def classify(ras,decs):
    """
    This is a function that applies selection criteria to sources with
    time series data in strip82 and identifies them as SNe, QSOs, variable stars,
    or standard stars using an sql query.

    Args:
        ras (2D array): An N row array containing ra values for the sources
        decs (2D array): An N row array containing dec values for the sources
    Returns:
        out_array(2D array): A N row array that tells you which objects are either
        0 = SN, 1 = QSO, 2 = variable star, and 3 = standard star.
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Clsssify astronomical sources with time series data in stripe82
    #
    #
    #
    # CALLING SEQUENCE:
    #   classify(ras,decs)
    #
    # INPUTS:
    #   ras - A 2D numpy array with shape (N,1).
    #   decs - A 2D numpy array with shape (N,1).
    #-------------------------------------------------------------

    # AB start timer
    time0 = time.time()
    #AB initalize arrays for loops
    utest = [] #u mag
    gtest = [] #g mag
    rtest = [] #r mag
    itest = [] #i mag
    ztest = [] #z mag
    mjdtest = [] #mjd
    errtest = [] #error in g mag

    #AB make 2D input array into 1D
    ras = np.hstack(ras)
    decs = np.hstack(decs)
    #AB convert to list for sql
    ras = ras.tolist()
    decs = decs.tolist()
    #AB loop for every source
    for i in range(len(ras)):
        #AB sql query that takes ra and dec of each source
        os.system("python ../../week13/azeem_stripe82query.py" + " " + str(ras[i]) + " " + str(decs[i]) + " " + ' > sqlresult' + str(i)+ '.data')
        #AB generate .dat file containing ra, dec, ugriz mags, g error, mjd, and distance columns
        ra,dec,u,gmag,r,imag,z,err,mjd,dist = np.genfromtxt('sqlresult' + str(i)+ '.data', unpack=True, skip_header=2)
        #AB make these into arrays to use
        utest.append(u)
        gtest.append(gmag)
        rtest.append(r)
        itest.append(imag)
        ztest.append(z)
        mjdtest.append(mjd)
        errtest.append(err)
    err = np.array(errtest)
    u = np.array(utest)
    gmag = np.array(gtest)
    r = np.array(rtest)
    imag = np.array(itest)
    z = np.array(ztest)
    mjd = np.array(mjdtest)

    #AB initalize output array with all 0s
    out_array = np.zeros(len(ras))
    #AB loop over every source to see which objects meet the criteria below
    for i in range(len(ras)):
        #AB don't want coadd
        mjd1 = mjd[i]
        good = np.where(mjd1>0)
        #AB want our variables to meet the condition above
        gmag1 = gmag[i][good]
        mjd1 = mjd[i][good]
        err1 = err[i][good]
        u1 = u[i][good]
        r1= r[i][good]
        i1= imag[i][good]
        z1= z[i][good]

        #AB don't want extraneous data points
        bad = np.where((err1<0.5) & (gmag1<np.median(gmag1)+1))
        g1 = gmag1[bad]
        mjd1 = mjd1[bad]
        err1 = err1[bad]
        u1 = u1[bad]
        r1= r1[bad]
        i1= i1[bad]
        z1= z1[bad]

        #AB now for our criteria
        #AB it is a SN if it has a standard deviation in g > 0.5
        if np.std(g1)>0.5:
            out_array[i] = 0
        #AB it is a QSO if it meets these color cuts from HW5 and has .15<std in g<0.5
        if u1[i]-r1[i]<1 and u1[i]-r1[i]>-0.09 and g1[i]-i1[i]<0.7 and g1[i]-i1[i]>-0.26 and g1[i]-r1[i]<0.68 and g1[i]-r1[i]>-0.12 and u1[i]-g1[i]>-0.15 and u1[i]-g1[i]<0.71 and np.std(g1)>.15 and np.std(g1)<.5: #QSO
            out_array[i] = 1
        #AB it is a standard star if it has a final standard error of the mean r band mag < 0.044
        elif np.std(r1)/np.sqrt(len(ras)) < 0.044:
            out_array[i] = 3
        #AB it is a variable star if it has 0.1<std in g < 0.5
        elif np.std(g1)>0.1 and np.std(g1)<0.5: #Var star
            out_array[i] = 2
    time1 = time.time()
    out_array = np.array(out_array)
    out_array = np.vstack(out_array)
    print("    the total number of seconds code takes to run [s]: {0:.3f}".format((time1-time0)))

    return out_array

#AB now feed in your 2D ra and dec arrays
if __name__ == '__main__':
    out_array = classify(ras,decs)
    print(out_array)
