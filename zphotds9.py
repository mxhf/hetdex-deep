
from astropy.io import fits
from astropy.wcs import WCS
import sys
from astropy.table import Table
import os

fncube = sys.argv[1]
hdu = fits.open(fncube)
fncatalog = "../pdz_cosmos2015_v1.3.fits"
fnout=fncube.replace(".fits.gz","_photz.reg")

if os.path.exists(fnout):
    print("{} exists already. Aborting.".format(fnout))
    sys.exit(1)

print("Reading {} ...".format(fncatalog))
t2 = Table.read( fncatalog )

wcs = WCS(hdu[0].header)

wcs = wcs.dropaxis(2)

ra_max,dec_min = wcs.wcs_pix2world(0,0,0)
ra_min,dec_max = wcs.wcs_pix2world(hdu[0].data.shape[1],hdu[0].data.shape[0],0)

jj  = (t2["RA"] > ra_min) * (t2["RA"] < ra_max)
jj *= (t2["DEC"] > dec_min) * (t2["DEC"] < dec_max)

sum(jj)

s = """\
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
"""

for r in t2[jj]:
    s += "point({:.6f},{:.6f}) # point=circle color=red font=\"helvetica 6 normal roman\" text={{z = {:.2f}}} \n".format(r["RA"],r["DEC"],r["Z_MED_PDZ"])

print("Writing {} ...".format(fnout))
with open(fnout,'w') as f:
    f.write(s)
    print("Wrote {}.".format(fnout))
