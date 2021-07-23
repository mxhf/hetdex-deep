# build catalog

from astropy.table import Table

import spectrum
from collections import OrderedDict
import argparse
from astropy.io import fits
from astropy import wcs
import argparse
import sys
import numpy as np
from scipy import interpolate

from astropy.io import ascii

from scipy.optimize import least_squares

import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

#generic gaussian
def gauss(sigma, xsq, norm=True):
    if norm:
        return 1./(sigma * np.sqrt(2. * np.pi) ) * np.exp( -xsq/(2. * sigma**2.))
    else:
        return np.exp( -xsq/(2. * sigma**2.))


def gauss3D(xc, yc, zc, sigmaxy, sigmaz, XX, YY, ZZ):
    RRSQ = (XX - xc)**2 + (YY - yc)**2 
    g = gauss(sigmaxy, RRSQ, norm=True) * gauss(sigmaz, (ZZ - zc)**2., norm=True)
    return g


fit_iteration_count = 0

def peval(p, XX, YY, ZZ):
    A, xc, yc, zc, sigmaxy, sigmaz = p
    g = gauss3D(xc, yc, zc, sigmaxy, sigmaz, XX, YY, ZZ) /  2.5066
    #g = g/np.sum(g)
    model = A * g
    #print("X", np.sum(g), np.sum(model))
    return model

def resid(p, XX, YY, ZZ, CC, EE, ii, DEBUG=False): 
    global fit_iteration_count
    fit_iteration_count += 1
    model = peval(p, XX[ii], YY[ii], ZZ[ii])
    A, xc, yc, zc, sigmaxy, sigmaz = p
    res = (CC[ii] - model)/EE[ii]
    if DEBUG:
        print(fit_iteration_count, A, xc, yc, zc, sigmaxy, sigmaz, np.sum(res**2.) /np.sum(ii) )
    return res
    
def rchisq(p, XX, YY, ZZ, CC, EE, ii, model=None):
    """
    Compute reduced chisq.
    """
    if type(model) == type(None):
        _model = peval(p, XX[ii], YY[ii], ZZ[ii])
    else:
        _model = model[ii]
    A, xc, yc, zc, sigmaxy, sigmaz = p
    res = (CC[ii] - _model)/EE[ii]
    return np.sum(res**2.)/(np.sum(ii) - len(p))


def fit_pointsource(id, x_com, y_com, z_com, xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, counts,\
                    CC, MM, EE, USEMAP, DEBUG, max_nfev=25):
    global fit_iteration_count

    if not USEMAP:
        print("Not using map.")
        ii  = XX >= xmin
        ii *= XX <= xmax
        ii *= YY >= ymin
        ii *= YY <= ymax
        ii *= ZZ >= zmin
        ii *= ZZ <= zmax
    else:
        ii = MM == id
        
    print("    Source {} has {} assigned spaxels.".format(id, np.sum(ii)))

    sigmaxy = np.sqrt(dx**2. + dy**2.)/1.414

    p0 = [counts, x_com, y_com, z_com, sigmaxy, dz]
    print("    Initial guesses are: ", p0)

    

    bestfit = least_squares(resid, p0, args=(XX, YY, ZZ, CC, EE, ii, DEBUG), max_nfev=max_nfev)


    model = peval(bestfit.x, XX, YY, ZZ)

    if DEBUG:
        # peak in pixel counts
        cnts_max = np.nanmax(CC[int(z_com)][ymin:ymax,xmin:xmax])
        N = zmax-zmin
        f = plt.figure(figsize=[20,10])
        for i, z in enumerate( range(zmin,zmax) ):
            plt.subplot(N, 3, i*3+1)
            plt.imshow( CC[z][ymin:ymax,xmin:xmax], vmin=0, vmax=cnts_max )

            plt.subplot(N, 3, i*3+2)
            plt.imshow( model[z][ymin:ymax,xmin:xmax], vmin=0, vmax=cnts_max )

            plt.subplot(N, 3, i*3+3)
            plt.imshow( CC[z][ymin:ymax,xmin:xmax] - model[z][ymin:ymax,xmin:xmax], \
                       vmin=-cnts_max/10., vmax=cnts_max/10.  )
        #plt.show()
        f.tight_layout()
        plt.savefig("fit_pointsource.pdf")
        
        print("bestfit", bestfit)
        
    return bestfit


#import os
#print("Limiting number of cores to 1.")
#os.environ['OPENBLAS_NUM_THREADS'] = '1' 
#os.environ['MKL_NUM_THREADS'] = '1'.

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--incube', type=str,                    
                    help='Input cube.')
parser.add_argument('-m', '--map', type=str,                    
                    help='Input map.')
parser.add_argument('-o', '--output_catalog', type=str, default='',
                    help='Output catalog.')
parser.add_argument('-c', '--calibration', type=str, default='',
                    help='Calibration.')
parser.add_argument('-C', '--fluxconversion', type=float, default=None,
                    help='Per count flux (alternative to --calibration) [erg/s/cm^2/A/cnt].')
parser.add_argument('-n', '--noisemodel', type=str, default='',                    
                    help='Detect noise model')

parser.add_argument('--dxmax', type=float, default = 4.,                  
                    help='X size threshold for point source fitting in px.')
parser.add_argument('--dymax', type=float, default = 4.,                  
                    help='Y size threshold for point source fitting in px.')
parser.add_argument('--dzmax', type=float, default = 10.,                  
                    help='Z size threshold for point source fitting in px.')

parser.add_argument('--usemap', action='store_true',                  
                    help='Use flag to only is spaxel that are flagged in' \
                         'the segmentation map rather than box aournd source.')

parser.add_argument('--id', type=float, default = None,                  
                    help='Only run for this id (for debugging).')

parser.add_argument('--debug', action='store_true',                  
                    help='Debug mode.')

parser.add_argument('--overwrite', action='store_true',                  
                    help='Overwrite existing output file.')

args = parser.parse_args(sys.argv[1:])

# Debug mode
DEBUG = args.debug
if DEBUG:
    from matplotlib import pyplot as plt
    
if not args.overwrite:
    if os.path.exists(args.output_catalog):
        print("Ouput file {} exists, aborting.".format(args.output_catalog))
        sys.exit(1)
    

# fit gaussians to pointsources
FIT_POINTSOURCES = True

names   = [ "id", "N", "counts", "flux", "ra_com", "dec_com", "ddec", "dra", "x_com", "y_com", "z_com", 
            "dx", "dy", "dz", "sqrt(ev1)", "sqrt(ev2)",
            "size2d",
            "x_ext", "y_ext", "z_ext",
            "wl_com", "dwl", "xmin","xmax", "ymin", "ymax", "zmin", "zmax",
            "psfit_cnts", "psfit_flux", "psfit_xc", "psfit_yc", "psfit_zc", \
            "psfit_ra", "psfit_dec", "psfit_wl",\
            "psfit_sigxy", "psfit_sigz", \
            "psfit_dxy", "psfit_dwl", \
            "psfit_rchisq",
            "rms_wl_com"
          ]

t = Table(names=names)

dtypes   =  {"id" : int, "N" : int, "counts" : float, "flux" : float, 
             "ra_com" : float, "dec_com" : float, 
             "dra" : float, "ddec" : float,
             "x_com": float, "y_com": float, "z_com": float,
             "dx": float, "dy": float, "dz": float, "sqrt(ev1)": float, "sqrt(ev2)": float,
             "size2d": float,
             "x_ext": int, "y_ext": int, "z_ext": int, 
             "wl_com" : float, "dwl": float, 
             "xmin": int,"xmax": int, "ymin": int, "ymax": int, "zmin": int, "zmax": int,
            "psfit_cnts" : float, "psfit_flux" : float, "psfit_xc" : float, "psfit_yc" : float, "psfit_zc" : float, 
             "psfit_ra" : float, "psfit_dec" : float, "psfit_wl" : float,
             "psfit_sigxy" : float, "psfit_sigz" : float, 
             "psfit_dxy" : float, "psfit_dwl" : float, 
             "psfit_rchisq" : float,
             "rms_wl_com" : float
            }


formats  =  {"id" : "5d", "N" : "5d", "counts" : ".4e", "flux" : ".4e", "ra_com" : "12.6f", "dec_com" : "12.6f", 
             "dra" : "3.1f", "ddec" : "3.1f", 
             "x_com": "8.2f", "y_com": "8.2f", "z_com": "8.2f", 
             "dx": "8.2f", "dy": "8.2f", "dz": "8.2f", "sqrt(ev1)": "8.2f","sqrt(ev2)": "8.2f",
             "size2d": "8.2f",
             "x_ext": "4d", "y_ext": "4d", "z_ext": "4d", 
             "wl_com" : "8.2f", "dwl": "8.2f", "xmin": "4d","xmax": "4d", "ymin": "4d", "ymax": "4d", "zmin": "4d", "zmax": "4d",
             "psfit_cnts" : ".4e", "psfit_flux" : ".4e", "psfit_xc" : "8.2f", "psfit_yc" : "8.2f", "psfit_zc" : "8.2f", 
             "psfit_ra" : "12.6f", "psfit_dec" : "12.6f", "psfit_wl" : "8.2f",
             "psfit_sigxy" : "8.2f", "psfit_sigz" : "8.2f", 
             "psfit_dxy" : "8.2f", "psfit_dwl" : "8.2f",
             "psfit_rchisq" : "8.2f",
             "rms_wl_com" : ".4e"}

units    =  {"id" : "", "N" : "px", "counts" : "counts", "flux" : "erg s^-1 cm^-2", 
             "ra_com" : "RA[J2000]", "dec_com" : "Dec[J2000]", "dra" : "arcsec", "ddec" : "arcsec", 
             "x_com": "px", "y_com": "px", "z_com": "px", 
             "dx": "px", "dy": "px", "dz": "px", "sqrt(ev1)": "px", "sqrt(ev2)": "px", 
             "size2d": "px",
             "x_ext": "px", "y_ext": "px", "z_ext": "px", 
             "wl_com" : "A", "dwl": "A", "xmin": "px","xmax": "px", "ymin": "px", "ymax": "px", "zmin": "px", "zmax": "px",
             "psfit_cnts" : "counts", "psfit_flux" : "erg s^-1 cm^-2", "psfit_xc" : "px", "psfit_yc" : "px", "psfit_zc" : "px", 
             "psfit_ra" : "RA[J2000]", "psfit_dec" : "Dec[J2000]", "psfit_wl" : "A",
             "psfit_sigxy" : "px", "psfit_sigz" : "px", 
             "psfit_dxy" : "arcsec", "psfit_dwl" : "A",
             "psfit_rchisq" : "",
             "rms_wl_com" : "counts"}

print("Building catalog for {} ...".format(args.incube))

print("Reading datacube {} ...".format(args.incube))
s = spectrum.readSpectrum(args.incube)
c = s.data
wlgrid = s.grid()
    

if FIT_POINTSOURCES:
    if args.noisemodel == "":
        print("You must specify a noise model for pointsource fitting.")
        sys.exit(1)
    print("Reading noisemodel {} ...".format(args.noisemodel))
    with open(args.noisemodel) as f:
        ll = f.readlines()
    pnoisemodel = [float(l) for l in ll]

    # build error cube
    # I am sure there is better ways of doing this...
    # here we evaluate the noise model at every wavelngth (= slice) and 
    # create a cube where each slice contans that value
    WW = np.ones_like(c)

    for i in range(len(wlgrid)):
        WW[i] = WW[i] * wlgrid[i]

    EE = np.polyval(pnoisemodel,WW)

print("Reading map {} ...".format(args.map))
outmap = fits.getdata(args.map)

    
zz,yy,xx = [np.arange(s, dtype=int) for s in c.shape]
YY,ZZ,XX = np.meshgrid(yy, zz, xx)

w = wcs.WCS(s.hdu)
platescale = s.hdu.header["CDELT2"]


if args.fluxconversion != None and args.calibration != "":
    print("ERROR: You can't specify a flux conversion file AND a constant count to flux conversion at the same time.")
    sys.exit(1)
    
if args.calibration != "":
    # load flux calibration
    print("Using thsi spectral response function {}".format(args.calibration))
    ftcal = args.calibration
    tcal = ascii.read(ftcal, format="fixed_width")
    cal_interp = interpolate.interp1d(tcal["wl[A]"], tcal["cal[erg/s/cm^2/A/cnt]"], fill_value='extrapolate')
elif args.fluxconversion != None:
    print("Using constant count to flux conversion of {} erg/s/cm^2/A/cnt".format(args.fluxconversion))
    cal_interp = interpolate.interp1d(wlgrid, [args.fluxconversion]*len(wlgrid) )
else:
    print("ERROR: You can specify either flux conversion file or a constant count to flux conversion.")

spec_response = cal_interp(wlgrid)

for n in names:
    t[n].dtype  = dtypes[n]
    t[n].format = formats[n]
    t[n].unit =  units[n]

 
c_dummy = np.zeros_like(c)

rr = np.sort( np.unique( outmap.flatten() ) )
print ("{} segments in total.".format(np.sum(rr>0)))

# Generate catalog
for r in rr[rr>0]:
    if args.id != None:
        if r != args.id:
            continue
            
    print("Source {} ...".format( r ))
        
    ii = outmap == r
    N = np.sum(ii)
    M = np.sum(c[ii])
    c_dummy[:,:,:] = 0.
    c_dummy[ii] = c[ii]
    
    
    # integrate flux
    MperWL = np.sum(np.sum(c_dummy, axis=2), axis=1)
    flux = ( np.sum(MperWL * spec_response) ) * s.step

    
    x_com = np.sum(c[ii]*XX[ii])/M
    y_com = np.sum(c[ii]*YY[ii])/M
    z_com = np.sum(c[ii]*ZZ[ii])/M
    
   
    mx = np.max(c_dummy)
    scale = 10000.
    
    #wl_com = float(ww_interp(z_com))
    zmax,zmin = (ZZ[ii].max()) , (ZZ[ii].min())
    ymax,ymin = (YY[ii].max()) , (YY[ii].min())
    xmax,xmin = (XX[ii].max()) , (XX[ii].min())
    
    #cov3d = np.cov( [XX.flatten(), YY.flatten(), ZZ.flatten()], fweights=np.array( c_dummy.flatten()*scale/mx , dtype=int) )
    #cov2d = np.cov( [XX.flatten(), YY.flatten()], fweights=np.array( c_dummy.flatten()*scale/mx , dtype=int) )
    #print("1 cov2d, cov3d = ", cov2d, cov3d)
    # faster version
    subdata = c_dummy[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten()
    fweights=np.array( subdata*scale/mx , dtype=int)
    if len(subdata) > 1 and any(fweights > 0):
        fweights[fweights<0] = 0
        cov3d = np.cov( [XX[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten(), \
                         YY[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten(),
                         ZZ[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten()],fweights=fweights\
                        )

        #cov2d = np.cov( [XX[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten(), YY[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten()],\
        #               fweights=np.array( c_dummy[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].flatten()*scale/mx , dtype=int) )
        cov2d = cov3d[:2,:2]
        #print("2 cov2d, cov3d = ", cov2d, cov3d)
        #1/0

        sqr_cov3d = np.sqrt( cov3d)
        dx = sqr_cov3d[0,0]
        dy = sqr_cov3d[1,1]
        dz = sqr_cov3d[2,2]
        #print(dx,dy,dz)

        eigenval,eigenvec = np.linalg.eig(cov2d)
        #print(np.sqrt( eigenval ))
        #print(np.sqrt(np.sum(eigenval)))
        e1,e2 = np.sqrt( eigenval ) 
    else:
        dx = 1.
        dy = 1.
        dz = 1.
        e1 = 1.
        e2 = 1.
        
        
    
    #dx = np.sqrt( np.sum( c[ii] * (XX[ii] - x_com)**2. ) / M ) * 2.35 # FWHM
    #dy = np.sqrt( np.sum( c[ii] * (YY[ii] - y_com)**2. ) / M ) * 2.35 # FWHM
    #dz = np.sqrt( np.sum( c[ii] * (ZZ[ii] - z_com)**2. ) / M ) * 2.35 # FWHM
    
    size2d = np.sqrt(np.sum(eigenval))  # FWHM

    z_ext = zmax - zmin
    x_ext = xmax - xmin
    y_ext = ymax - ymin
    
    # use atropy wcs to convert x,y,z, dz to ra,dec,wl  ... and in bit of a hack ... dwl
    voxel = w.wcs_pix2world([[x_com,y_com,z_com],[x_com,y_com,z_com+dz] ],1)
    dwl = np.abs(voxel[1][2] - voxel[0][2])
    ra_com, dec_com, wl_com = voxel[0]
    
    dra = platescale*dx * 3600.
    ddec = platescale*dy * 3600.
    
    # default nan for PSF fits values if fit is not carried out
    psfit_cnts, psfit_flux, psfit_xc, psfit_yc, psfit_zc = np.nan, np.nan, np.nan, np.nan, np.nan
    psfit_ra, psfit_dec, psfit_wl = np.nan, np.nan, np.nan
    psfit_sigxy, psfit_sigz = np.nan, np.nan
    psfit_dxy, psfit_dwl = np.nan, np.nan
    psfit_rchisq  = np.nan
    
    rms_wl_com = np.polyval(pnoisemodel,wl_com)
    
    if FIT_POINTSOURCES:
        print("Fitting gaussian model to source {} ...".format(r))
        
        if (xmax-xmin) < 4:
            print("   x extend is less than 5. Skipping ....")
        elif (ymax-ymin) < 4:
            print("   y extend is less than 5. Skipping ....")
        elif (zmax-zmin) < 4:
            print("   z extend is less than 5. Skipping ....")
        elif dx > args.dxmax:
            print("   dx = {} is larger than threshold {}. Skipping ....".format(dx, args.dxmax))
        elif dy > args.dymax:
            print("   dy = {} is larger than threshold {}. Skipping ....".format(dy, args.dymax))
        elif dz > args.dzmax:
            print("   dz = {} is larger than threshold {}. Skipping ....".format(dz, args.dzmax))
        else:    

            # run the pointsource fit
            #if True:
            try:
                #print("fit_pointsource",r, x_com, y_com, z_com, xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, M,\
                #            c, outmap, EE, args.usemap, DEBUG)
                bestfit = fit_pointsource(r, x_com, y_com, z_com, xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, M,\
                            c, outmap, EE, args.usemap, DEBUG)
                
                print("    bestfit.success :", bestfit.success)
                
                if not bestfit.success:
                    print("    Fit failed.")
                else:
                    p = bestfit.x
                    # print("P: ", p)
                    psfit_cnts, psfit_xc, psfit_yc, psfit_zc, psfit_sigxy, psfit_sigz = p
                    psfit_sigxy = np.abs(psfit_sigxy) # make positive definite
                    psfit_sigz = np.abs(psfit_sigz) # make positive definite
                    # integrate model
                    model = peval(p, XX, YY, ZZ)
                    MperWL = np.sum(np.sum(model, axis=2), axis=1)
                    psfit_flux =  ( np.sum( MperWL * spec_response) ) * s.step

                    psfit_rchisq = rchisq(p, XX, YY, ZZ, c, EE, ii, model=model)

                    print("   Bestfit parameters after {} iterations are: ".format(fit_iteration_count), p)
                    print("   Reduced chisq is: {:.4f}".format(psfit_rchisq) )

                    psfit_dxy = platescale * psfit_sigxy * 3600.


                    voxel = w.wcs_pix2world([[psfit_xc,psfit_yc,psfit_zc],[psfit_xc,psfit_yc,psfit_zc + psfit_sigz] ],1)
                    psfit_dwl = voxel[0][2] - voxel[1][2]
                    psfit_ra, psfit_dec, psfit_wl = voxel[0]
                    
                fit_iteration_count = 0
            except:
                print("    Fit failed.")

    
    t.add_row([ r, N, M, flux, ra_com, dec_com, dra, ddec, x_com, y_com, z_com, \
               dx, dy, dz, e1,e2, size2d, x_ext, y_ext, z_ext, wl_com, dwl, \
               xmin,xmax, ymin, ymax, zmin, zmax, \
               psfit_cnts, psfit_flux, psfit_xc, psfit_yc, psfit_zc, \
               psfit_ra, psfit_dec, psfit_wl,\
               psfit_sigxy, psfit_sigz, \
               psfit_dxy, psfit_dwl, \
               psfit_rchisq,\
               rms_wl_com] )
               
    
    
t.write(args.output_catalog, format="ascii.ecsv", overwrite=True)

print("Wrote {}.".format(args.output_catalog))

