#!/usr/bin/env python
import os
import numpy as np
import sys
import argparse
from collections import OrderedDict

from astropy.table import Table, Column
from astropy.stats import biweight_location
from astropy.io import fits
from astropy.io import ascii

from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve

import spectrum


# simple continuum removal though spline interpolation
def confitSpl(wls, s, n = 10, kappaL=2.5, kappaU=2.5, output_fit=False, smooth=0., mask=None, PLOT=False, maxiter=15):
    if mask == None:
        mask_ = wls > -1 # use all
    else:
        mask_ = mask.copy()
    mask_ *= ~ np.isnan(s)
    l = len(wls)
    l = np.floor(l/n)
    dwl = (wls[-1]-wls[0])/n

    niter = 0
    nmasked = len(mask_[~mask_])
    while niter < maxiter:
        bwls = []
        bs   = []

        # put one point at the blue end, window only half the normal binsize
        wlstart = wls[0]
        wlstop  = wls[0] + dwl/2.
        ii = (wls >= wlstart) * (wls <= wlstop)
        if type(mask_) != type(None): ii *= mask_

        binned_wls = np.mean( wls[ii] )
        bwls.append( binned_wls )
        bs.append(    np.median(   s[ii] ) )
        # normal points, normal binsize
        for i in range(n-1):
                wlstart = wls[0]  + dwl/2. + dwl * i
                wlstop  = wls[0]  + dwl/2. + dwl * (i + 1)
                ii = (wls >= wlstart) * (wls <= wlstop)
                if type(mask_) != type(None): ii *= mask_
                binned_wls = np.mean( wls[ii] )
                bwls.append( binned_wls )
                bs.append(    np.median(   s[ii] ) )
        # put one point at the red end, window only half the normal binsize
        wlstart = wls[-1] - dwl/2.
        wlstop  = wls[-1]
        ii = (wls >= wlstart) * (wls <= wlstop)
        if type(mask_) != type(None): ii *= mask_
        binned_wls = np.mean( wls[ii] )
        bwls.append( binned_wls )
        bs.append(    np.median(   s[ii] ) )

        bwls = np.array( bwls )
        bs = np.array( bs )
        jj = ~np.isnan(bwls) * ~np.isnan(bs)
        tck = interpolate.splrep(bwls[jj],bs[jj],s=smooth)
        c = interpolate.splev(wls,tck,der=0)

        res = s-c
        sigma = np.std(res[mask_])

        inliers  = ( res) <= kappaU*sigma
        inliers *= (-res) <= kappaL*sigma

        mask_ *= inliers
        nmasked_new = len(mask_[~mask_])
        if nmasked_new == nmasked:
            break
        nmasked = nmasked_new

        niter += 1
    if PLOT:
        f=plt.figure()
        plt.plot(wls,s)
        plt.plot(wls,c)
        plt.plot(wls[~mask_],s[~mask_],'r.')
        plt.ylim([-1.,1.])

    # filter lowest and highest 3 fourier channels
    sc = s-c

    if output_fit:
        return sc,c
    else:
        return sc 


def line_detect(ww, csout, threshold):
    # line detection (everything above cetain threshold)
    jj = csout > threshold

    # labelling line detections
    label = skimage.morphology.label(jj)
    ll = np.unique( label )

    lineset = []
    dlineset = []

    for l in ll:
        if l == 0:
            continue
        ii = l == label
        f = np.sum( csout[ii] )
        wl_com = np.sum( ww[ii]*csout[ii] ) /np.sum(csout[ii] )
        #print("{} {:.2f}A {:.2f}".format(l, wl_com, f))
        lineset.append(wl_com)
        dlineset.append(2.)
    return lineset, dlineset, jj


def masked_biweight(cube_slice, mask):
    return biweight_location( cube_slice[mask] )


def masked_sum(cube_slice, mask):
    return np.sum( cube_slice[mask] )


def extract(r, s, outmap, method=masked_biweight):
    mask = np.sum( outmap == r['id'], axis=0) > 0

    sout = np.zeros( s.data.shape[0]  )
    N = np.sum(mask)
    for i in range(s.data.shape[0]):
        #sout[i] = np.sum( s.data[i][mask] )
        #sout[i] = biweight_location( s.data[i][mask] )*N # statistically stable mean times number of pixel
        sout[i] = method( s.data[i], mask )

    ww = s.grid()
    return ww,sout, mask


def nextract(r, s, outmap, method=masked_biweight, MAX_SAMPLES = 30):
    """
    Shifts the aperture corresponding to a detection around in the noise cube to
    sample the noise N times. The mask is shifted such as to not overlap
    wiht he previous location of the mask

    Original, actual mask
    xxxxxxxxxx
    xxx..xxxxx
    xx...xxxxx
    xxx..xxxxx
    xxxxxxxxxx
    xxxxxxxxxx

    Sample 1
    x..xxxxxxx
    ...xxxxxxx
    x..xxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx

    Sample 2
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx
    x..xxxxxxx
    ...xxxxxxx
    x..xxxxxxx

    Sample 3
    xxxx..xxxx
    xxx...xxxx
    xxxx..xxxx
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx

    Sample 4
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx
    xxxx..xxxx
    xxx...xxxx
    xxxx..xxxx

    Sample 5
    xxxxxxx..x
    xxxxxx...x
    xxxxxxx..x
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx

    Sample 6
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxxxxx
    xxxxxxx..x
    xxxxxx...x
    xxxxxxx..x
    """

    mask = np.sum( outmap == r['id'], axis=0) > 0

    # determine masked, region size
    xx = np.arange(mask.shape[1])
    yy = np.arange(mask.shape[0])

    X,Y = np.meshgrid(xx,yy)
    minx,maxx = X[mask].min(), X[mask].max()
    miny,maxy = Y[mask].min(), Y[mask].max()
    sx = maxx - minx + 2
    sy = maxy - miny + 2
    nx,ny = mask.shape[1]//sx, mask.shape[0]//sy

    #f = plt.figure(figsize = [7,7])
    #plt.imshow(mask)

    all_sout = []
    count = 1
    print("Sampling noise in {} separate regions".format( min(nx * ny, MAX_SAMPLES)  ) )
    for i in range(nx):
        if count > MAX_SAMPLES:
            break
        
        for j in range(ny):
            print("Sampling noise column {} and row {}".format(i, j))
            _mask = mask < 0 # empty mask
            #print(sy*(j), sy*(j+1),sx*(i), sx*(i+1), miny, miny+sy, minx, minx+sx )

            m = mask[miny:miny+sy, minx:minx+sx]
            _mask[sy*(j):sy*(j)+m.shape[0],sx*(i):sx*(i)+m.shape[1]] = m
            #f = plt.figure(figsize = [7,7])
            #plt.imshow(_mask)

            sout = np.zeros( s.data.shape[0]  )
            for k in range(s.data.shape[0]):
                sout[k] = method( s.data[k], _mask )
            all_sout.append(sout)
            count += 1

    ww = s.grid()
    return ww, np.array(all_sout), mask


# In[32]:


def save_spec(ifu, field, id, sout, nsout=None, mask=None, DEBUG=False):
    t = Table()
    for dt in ['D', 'D', '16A', 'D', 'D', 'D', 'D', 'D', 'D']:
        t.add_column(Column(dtype='d'))

    t.add_row()

    hdu1 = PrimaryHDU()
    hdu2 = BinTableHDU(t)

    hdu3 = ImageHDU()
    hdu4 = ImageHDU()
    hdu5 = ImageHDU()

    hdu3.data = [sout]
    
    if type(nsout) == type(None):
        # if no noise spectrum is given save just array of nans
        hdu4.data = [[np.nan]*len(sout)]
    else:    
        mm = np.std(nsout,axis=0)
        hdu4.data = [mm]

    if type(mask) == type(None):
        # if no maks is given, mask nan values in spectrum and noise spectrum
        hdu5.data = np.array( [~np.isnan(mm) * ~np.isnan(sout)], dtype=int)
    else:
        hdu5.data = mask
    
    new_hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])
    
    spec_filename = '../data/extractions/{}/{}_{:05d}.fits'.format(field, ifu, id)
    h,_ = os.path.split(spec_filename)
    if not os.path.exists(h):
        os.makedirs(h)
        
    new_hdul.writeto(spec_filename, overwrite=True)
    
    return spec_filename


def extract_all(r, c,  outmap, sraw, ns, MAX_SAMPLES, DEBUG=False):
    ww, sout, mask    = extract(r, s, outmap)
    if not type(sraw) == type(None):
        ww, sout_unsmoothed, __  = extract(r, sraw, outmap, method=masked_sum)
    else:
        sout_unsmoothed = None
        
    if not type(ns) == type(None):
        ww, nsout, mask = nextract(r, ns, outmap, MAX_SAMPLES = MAX_SAMPLES)
    else:
        nsout = None
        
    return ww, sout, sout_unsmoothed, nsout


def fit_continuum(ww, sout, sout_unsmoothed, DEBUG=False):

    # continuum removal

    csout, cont            = confitSpl(ww, sout, output_fit=True, PLOT=False)
    
    if not type(sout_unsmoothed) == type(None):
        ii = ~np.isnan(sout_unsmoothed)
        csout_unsmoothed, cont_unsmoothed = confitSpl(ww, sout_unsmoothed, output_fit=True, PLOT=False)
    else:
        cont_unsmoothed = None

    return cont, cont_unsmoothed


def compute_cont_stats(ww, wl_com, cal_cont, Nwin = 3, DEBUG=False):

    dwl = ww[1] - ww[0]

    int_cont = np.nansum(cal_cont) * dwl
    if DEBUG:
        print("Integrated continuum flux is {:.2e}".format(int_cont))
    mn_cont = np.nanmean(cal_cont) 
    if DEBUG:
        print("Mean continuum flux is {:.2e}".format(mn_cont))

    dd = np.abs(ww - wl_com)
    imin = np.argmin(dd)

    loc_cont = (cal_cont)[imin]
    if DEBUG:
        print("Continuum flux at detection wavelength {:.1f}A is {:.2e}".format(wl_com, loc_cont))

    win_cont = []

    for i in range(Nwin):
        iwlmin = int( len(ww) * (i)/Nwin)
        iwlmax = min(int( len(ww) * (i+1)/Nwin), len(ww)-1 )
        wlmin, wlmax = ww[iwlmin], ww[iwlmax]
        
    
        ii = (ww >= wlmin) * (ww < wlmax)
        mn_cont_win = np.nanmean((cal_cont)[ii]) 
        if DEBUG:
            print("Mean continuum flux in wl window {:.1f}A - {:.1f}A is {:.2e}".format(wlmin, wlmax, mn_cont_win))
        win_cont.append(mn_cont_win)

    return int_cont, mn_cont, loc_cont, win_cont


def prepend(path, prefix):
    h,t = os.path.split(path)
    return os.path.join(h, prefix + t)

    
def main(s, outmap, t, outcatalog, ids = None, ns=None, sraw=None, Nwin=3, MAX_SAMPLES=30, DEBUG=False):
    # definition of additional catalog columns to hold continuum information
    cc = OrderedDict()
    cc[ "int_cont" ] = { 'dtype' : float, 'unit' : "erg s^-1 cm^-2", 'fmt' : ".4e"} 
    cc[ "mn_cont" ]  = { 'dtype' : float, 'unit' : "erg s^-1 cm^-2 A^-1", 'fmt' : ".4e"} 
    cc[ "loc_cont" ] = { 'dtype' : float, 'unit' : "erg s^-1 cm^-2 A^-1", 'fmt' : ".4e"} 
    for i in range(Nwin):
        cc[ "win_cont_{}".format(i+1) ] = { 'dtype' : float, 'unit' : "erg s^-1 cm^-2 A^-1", 'fmt' : ".4e"} 

    # ann new columns to table
    for c in cc:
        col = Column(name=c, dtype=cc[c]['dtype'], unit=cc[c]['unit'], format=cc[c]['fmt'], data=[np.nan]*len(t))
        t.add_column(col)
            
    for id in ids:
        try:
            ii = t['id'] == id
            if DEBUG:
                print("Id {}".format(id))
            r = t[ii][0]
            wldetect = r["wl_com"]


            id, ra,dec = r["id"], r["ra_com"], r["dec_com"]

            # spectral extraction
            ww, sout, sout_unsmoothed, nsout = extract_all(r, c,  outmap, sraw, ns, MAX_SAMPLES, DEBUG=DEBUG)

            # fit continuum
            cont, cont_unsmoothed = fit_continuum(ww, sout, sout_unsmoothed, DEBUG=DEBUG)

            pixelscale = np.abs( s.hdu.header["CDELT1"]*3600. )
            A = pixelscale**2. # pixelsize

            #c = s.data
            if DEBUG:
                print("Using constant count to flux conversion of {} erg/s/cm^2/A/cnt".format(1e-17))
            cal_interp = interpolate.interp1d(ww, [1e-17]*len(ww) )

            # flux calibrated spectra
            cal_spec = sout*cal_interp(ww)*1e17/A
            #if not type(nsout) == type(None):
            #    cal_noisespec = nsout*cal_interp(ww)*1e17/A
            cal_cont = cont*cal_interp(ww)/A

            # compute statistics on continuum
            int_cont, mn_cont, loc_cont, win_cont = compute_cont_stats(ww, r["wl_com"], cal_cont, Nwin, DEBUG=DEBUG)

            t["int_cont"][ii] = int_cont
            t["mn_cont"][ii]  = mn_cont
            t["loc_cont"][ii] = loc_cont
            for i in range(Nwin):
                t[ "win_cont_{}".format(i+1) ] = win_cont[i]
        except:
            print("Continuum level estimation failed for source {}.".format(id))
    
    #outfilename = prepend( fcatalog , "cnt")
    t.write( outcatalog, format="ascii.ecsv", overwrite=True)
    print("Wrote {}.".format(outcatalog))
    
    #save_spec(ifu, field, id, sout, nsout, mm)


# In[30]:


if False:
    print("Start")

    IFU = "015"
    field = "gama09Efin"

    fcube = "../data/sf{}_{}.fits.gz".format(field, IFU)
    fcatalog = "../data/sf{}_{}.cat".format(field,IFU)

    # test
    fmap = "../data/map{}_{}.fits.gz".format(field, IFU)
    frawcube = "../data/{}_{}.fits.gz".format(field, IFU)
    fnoisecube = "../data/sfnc{}_{}.fits.gz".format(field, IFU)
    fcatalog = "../data/sf{}_{}.cat".format(field,IFU)


if True:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, metavar='infile',
                        help='Input cube.')
    parser.add_argument('--id', type=int, metavar='id',
                        help='Specific ID to work on (mostly for debugging).', default=None)
    parser.add_argument('-c', '--catalog', type=str, metavar='catalog',
                        help='Input catalog.')
    parser.add_argument('-m', '--map', type=str, metavar='map',
                        help='Input segmentation map.')
    parser.add_argument('-o', '--outcatalog', type=str, metavar='outcatalog',
                        help='Ouput catalog.')
    #parser.add_argument('--rawcube', type=str, metavar='rawcube',
    #                    help='Unsmoothed input cube.', default=None)
    args = parser.parse_args(sys.argv[1:])
    
    

fcube = args.infile
fmap = args.map
fcatalog = args.catalog
outcatalog = args.outcatalog

frawcube = None
fnoisecube = None
ids = None

# read tables
print("Reading catalog ... ")
t = ascii.read(fcatalog, format="ecsv")
    
# read spectra
print("Reading cubes ... ")
s = spectrum.readSpectrum(fcube)
outmap = fits.getdata(fmap)
if not frawcube == None:
    sraw = spectrum.readSpectrum(frawcube)
else:
    sraw = None
if not fnoisecube == None:
    ns = spectrum.readSpectrum(fnoisecube)
else:
    ns = None

# read tables
print("Reading catalog ... ")
t = ascii.read(fcatalog, format="ecsv")
  
if type(args.id) == type(None):
    ids = None
else:
    ids = [args.id]

if type(ids) == type(None):
    ids= np.sort(np.unique(t['id']))
    
if os.path.exists(outcatalog):
    print("{} already exists. Returning.".format(outcatalog))
    sys.exit(0)
 
main(s, outmap, t, outcatalog, ids = ids,      ns=ns, sraw=sraw, Nwin=3, MAX_SAMPLES=30, DEBUG=True)


