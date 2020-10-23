#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import spectrum
from astropy.io import ascii
from astropy.io import fits
from scipy import interpolate
import numpy
from matplotlib import pyplot as plot
import skimage.morphology
from scipy.optimize import least_squares
from astropy.stats import biweight_location
import numpy as np
from astropy.table import Table, Column


import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

import matplotlib.pyplot as plt

from astropy import wcs
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
    
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.units as u
import os
import pickle

from scipy.interpolate import interp1d


def register_ds9staircase():
    # register color map
    from matplotlib.cm import register_cmap, cmap_d

    colors = []
    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*.3,kk*.3,kk*1)  )

    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*.3,kk*1,kk*.3)  )
    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*1,kk*.3,kk*.3)  )
    colors = np.array(colors)
    xx = np.arange(len(colors), dtype=float)
    xx = xx/xx.max()

    ds9staircase = {'red': lambda v : np.interp(v, xx, colors[:,0]),
               'green': lambda v : np.interp(v, xx, colors[:,1]),
               'blue': lambda v : np.interp(v, xx, colors[:,2])}


    # Register all other colormaps
    register_cmap('ds9staircase', data=ds9staircase)
    
register_ds9staircase()




def show_collased(ax, s, r, cal_interp, width = 20./3600., vmin=0., vmax=8., colorbar=True):
    A = 0.5**2.
    print("show_slice new 2")
    ra,dec = r["ra_com"], r["dec_com"]

    hdu = s.hdu
    w = WCS(hdu.header)
    w = w.dropaxis(2)
    colormap = plt.get_cmap('ds9staircase')
    
    ww = s.grid()
    dwl = ww[1] - ww[0]
    sl = np.sum( hdu.data[r["zmin"]:r["zmax"]], axis=0)

    sl = sl * cal_interp( r["wl_com"]) * dwl/A * 1e18
    
    cos_term = np.cos(np.deg2rad(dec))

    sl[sl == 0] = np.nan
    #current_cmap = plt.cm.get_cmap()
    colormap.set_bad(color='grey')
    cax = ax.imshow(sl, vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest', cmap=colormap)
    x0,y0 = w.wcs_world2pix(ra - width/2./cos_term ,dec-width/2.,0)
    x1,y1 = w.wcs_world2pix(ra + width/2./cos_term, dec+width/2.,0)
    x,y = w.wcs_world2pix(ra,dec,0)
    ax.plot([x],[y],'x',c='white')
    
    if colorbar:
        cb = plt.colorbar(cax)
        cb.set_label('10$^{-18}$erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$')

        ax.set_xlim([x1,x0])
        ax.set_ylim([y0,y1])

        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    
    return sl

  

import glob
from astropy.io import fits
from astropy import wcs
import sys, traceback
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--field", type=str, default='COSMOSA',
                    help="Field to generate plots for.")
parser.add_argument("-i", "--ifu", type=str, default='022',
                    help="IFU to generate plots for.")
parser.add_argument("-j", "--object_id", type=float, default=np.nan,
                    help="If of object to plot.")

args = parser.parse_args(sys.argv[1:])


print("Start ...")

IFU = args.ifu


field = args.field
IFU = args.ifu
fincube = "data/sf2outcube_COSMOSA_{}.fits.gz".format(IFU)
foutmap = "data/map_{}_{}.fits.gz".format(field, IFU)
fincube = "data/sf2outcube_{}_{}.fits.gz".format(field, IFU)
frawincube = "data/outcube_{}_{}.fits.gz".format(field, IFU)
fnoisecube = "data/sf2ncoutcube_{}_{}.fits.gz".format(field, IFU)
ftcal = "mean_cal.txt"
object_id = args.object_id

# load flux calibration
tcal = ascii.read(ftcal, format="fixed_width") 
cal_interp = interp1d(tcal["wl[A]"], tcal["cal[erg/s/cm^2/A/cnt]"], kind='cubic', bounds_error=False)

fcatalog = "data/sf2outcube_COSMOSA_{}.cat".format(IFU)
t = ascii.read(fcatalog, format="ecsv")

ii = t['id'] == object_id
r = t[ii][0]

s = spectrum.readSpectrum(fincube)
wcs = WCS(s.hdu.header)
wcs = wcs.dropaxis(2)

f = plt.figure(figsize=[7,5])
ax = plt.subplot(111, projection=wcs)
ax.set_facecolor('grey')

width = 20./3600.
show_collased(ax, s, r, cal_interp, width=width)


plt.savefig("test.pdf")

sys.exit(0)


