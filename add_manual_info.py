# 
# does really nothing more than augmenting the atuomaticaly generated catalog (by build_catalog.py)
# and adds a comment column where the classification string from the manual catalog goes.

import sys
import numpy as np
import argparse
import os
from astropy.io import fits, ascii
import spectrum

from astropy.table import Column
#from astropy.wcs import WCS

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--map', type=str,                    
                    help='Input map.')
parser.add_argument('-c', '--manual_catalog', type=str, default='',
                    help='Manually generated catalog.')
parser.add_argument('-a', '--automated_catalog', type=str, default='',
                    help='Automaticaly generated catalog.')

args = parser.parse_args(sys.argv[1:])

catfile = args.manual_catalog 


# read manual catalog
with open(catfile, 'r') as f:
    ll = f.readlines()
    c = []
    for l in ll:
        if l.strip().startswith("#"):
            continue
        if l.strip() == "":
            continue
        tt = l.split()
        x,y,z = [int(np.round(float(f))) for f in tt[:3]]
        comment = ""
        if len(tt) > 3:
            #print(tt[3:], comment)
            comment = ','.join( tt[3:] )
        c.append([x,y,z, comment])
        
        
t = ascii.read(args.automated_catalog, format="ecsv")
if 'class' in t.colnames:
    t.remove_column('class')
if 'manualx' in t.colnames:
    t.remove_column('manualx')
if 'manualy' in t.colnames:
    t.remove_column('manualy')
if 'manualz' in t.colnames:
    t.remove_column('manualz')
t.add_column(Column(["NA"]*len(t), name='class', dtype='S100'))

t.add_column(Column([np.nan]*len(t), name='manualx', dtype=float))
t.add_column(Column([np.nan]*len(t), name='manualy', dtype=float))
t.add_column(Column([np.nan]*len(t), name='manualz', dtype=float))

m = spectrum.readSpectrum(args.map)
for r in c:
    ix = int( np.round(r[0]) -1. )
    iy = int( np.round(r[1]) -1.)
    iz = int( np.round(r[2]) -1.)
    segment_id = m.data[iz,iy,ix]
    ii = t['id'] == segment_id
    
    if segment_id < 1:
        print("Segment ID for manually flaggged source at {},{},{} is "\
              "0 i.e. has no corresponding segment in map.".format(ix,iy,iz))
        continue
        
    if sum(ii) != 1:
        print("Error, unable to find segment id {} in {}.".format(segment_id, args.automated_catalog))
        continue
        
    t['class'][ ii ] = r[3]
    
    t['manualx'][ ii ] = r[0]-1. 
    t['manualy'][ ii ] = r[1]-1. 
    t['manualz'][ ii ] = r[2]-1. 
    #print(r[3])

h,tail = os.path.split(args.automated_catalog)
#t.write( os.path.join(h,"a" + tail), format="ascii.ecsv",overwrite=True)
# overwrite input
t.write( args.automated_catalog, format="ascii.ecsv",overwrite=True)
print("Wrote {}".format(args.automated_catalog))
    
