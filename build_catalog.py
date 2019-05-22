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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--incube', type=str,                    
                    help='Input cube.')
parser.add_argument('-m', '--map', type=str,                    
                    help='Input map.')
parser.add_argument('-o', '--output_catalog', type=str, default='',
                    help='Output catalog.')
parser.add_argument('-c', '--calibration', type=str, default='',
                    help='Calibration.')

args = parser.parse_args(sys.argv[1:])


names   = [ "id", "N", "counts", "flux", "ra_com", "dec_com", "ddec", "dra", "x_com", "y_com", "z_com", 
           "dx", "dy", "dz", "sqrt(ev1)", "sqrt(ev2)",
           "size2d",
           "x_ext", "y_ext", "z_ext",
           "wl_com", "dwl", "xmin","xmax", "ymin", "ymax", "zmin", "zmax"]
t = Table(names=names)

dtypes   =  {"id" : int, "N" : int, "counts" : float, "flux" : float, 
             "ra_com" : float, "dec_com" : float, 
             "dra" : float, "ddec" : float,
             "x_com": float, "y_com": float, "z_com": float,
             "dx": float, "dy": float, "dz": float, "sqrt(ev1)": float, "sqrt(ev2)": float,
             "size2d": float,
             "x_ext": int, "y_ext": int, "z_ext": int, 
             "wl_com" : float, "dwl": float, 
             "xmin": int,"xmax": int, "ymin": int, "ymax": int, "zmin": int, "zmax": int}

formats  =  {"id" : "5d", "N" : "5d", "counts" : ".4e", "flux" : ".4e", "ra_com" : "12.6f", "dec_com" : "12.6f", 
             "dra" : "3.1f", "ddec" : "3.1f", 
             "x_com": "8.2f", "y_com": "8.2f", "z_com": "8.2f", 
             "dx": "8.2f", "dy": "8.2f", "dz": "8.2f", "sqrt(ev1)": "8.2f","sqrt(ev2)": "8.2f",
             "size2d": "8.2f",
             "x_ext": "4d", "y_ext": "4d", "z_ext": "4d", 
             "wl_com" : "8.2f", "dwl": "8.2f", "xmin": "4d","xmax": "4d", "ymin": "4d", "ymax": "4d", "zmin": "4d", "zmax": "4d"}

units    =  {"id" : "", "N" : "px", "counts" : "counts", "flux" : "erg s^-1 cm^-2", 
             "ra_com" : "RA[J2000]", "dec_com" : "Dec[J2000]", "dra" : "arcsec", "ddec" : "arcsec", 
             "x_com": "px", "y_com": "px", "z_com": "px", 
             "dx": "px", "dy": "px", "dz": "px", "sqrt(ev1)": "px", "sqrt(ev2)": "px", 
             "size2d": "px",
             "x_ext": "px", "y_ext": "px", "z_ext": "px", 
             "wl_com" : "A", "dwl": "A", "xmin": "px","xmax": "px", "ymin": "px", "ymax": "px", "zmin": "px", "zmax": "px"}

print("Building catalog for {} ...".format(args.incube))

s = spectrum.readSpectrum(args.incube)
c = s.data
wlgrid = s.grid()

zz,yy,xx = [np.arange(s, dtype=int) for s in c.shape]
YY,ZZ,XX = np.meshgrid(yy, zz, xx)

w = wcs.WCS(s.hdu)
platescale = s.hdu.header["CDELT2"]
outmap = fits.getdata(args.map)


# load flux calibration
ftcal = args.calibration
tcal = ascii.read(ftcal, format="fixed_width")
cal_interp = interpolate.interp1d(tcal["wl[A]"], tcal["cal[erg/s/cm^2/A/cnt]"], fill_value='extrapolate')


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
    voxel = w.wcs_pix2world([[x_com,y_com,z_com],[x_com,y_com,ZZ[ii].max()],[x_com,y_com,ZZ[ii].min()]],1)
    dwl = voxel[1][2] - voxel[2][2]
    ra_com, dec_com,wl_com = voxel[0]
    
    dra = platescale*dx * 3600.
    ddec = platescale*dy * 3600.
    t.add_row([ r, N, M, flux, ra_com, dec_com, dra, ddec, x_com, y_com, z_com, dx, dy, dz, e1,e2, size2d, x_ext, y_ext, z_ext, wl_com, dwl, xmin,xmax, ymin, ymax, zmin, zmax] )
    
t.write(args.output_catalog, format="ascii.ecsv", overwrite=True)

print("Wrote {}.".format(args.output_catalog))

