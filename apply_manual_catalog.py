#
# Takes manually build catalog and adds segments to segmentation map.

# coding: utf-8
import sys
import numpy as np
import argparse

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.wcs import WCS

import cloud_finder3 as cf

def apply_manual_catalog(c, m, d, wcs, p, sigma_grow_threshold, auto_and_man_label_start=5000, only_man_label_start = 10000):
    print("Maxlabel is: ", only_man_label_start)
    currentlabel = only_man_label_start

    count_labeled = 0
    s = "{:3s} {:3s} {:3s} {:3s} {:8s} {:7s} {:5s} {:6s} {:5s} {:6s}"
    print( s.format( "cid", "x", "y", "z", "val", "wl[A]", "val/noise", "N_labeled", "newid", "N_spax"))
    s = "==================================="
    print( s )

    new_regs = []

    
    
    for i,r in enumerate(c):
        x,y,z = r
        id = m[z-1,y-1,x-1]
        if id > 0:
            # was labeled before
            count_labeled += 1
            if  id < auto_and_man_label_start:
                m[m == id] += auto_and_man_label_start

        # was not labeled before
        else:
            print("Max Segment ID is: ", np.max(m))
            
            ra,dec,wl = wcs.wcs_pix2world(x,y,z,1)
            noise_at_wl = np.polyval(p,wl)
            val = d[z-1,y-1,x-1]
        
            threshold = sigma_grow_threshold * noise_at_wl
            print("threshold: ", threshold)

            # make sure als lest the manually labeled spaxel get a 
            # (if must signle-spaxel) segment
            m[z-1,y-1,x-1] += currentlabel+1
            # try growing the segment
            m = cf.grow_segment(d, m, x, y, z, threshold, currentlabel+1)

            N_labeled_new = np.sum(m > 0.) 
            print("currentlabel ", currentlabel)
            currentlabel += 1

            s = "{:3d} {:3d} {:3d} {:3d} {:6.3f} {:7.2f} {:5.1f} {:6d} {:5d} {:6d}"
            N_spax = np.sum(m == currentlabel)
            print( s.format( i, x,y,z,val,float(wl),val/noise_at_wl, N_labeled_new, currentlabel, N_spax)  )

            new_regs.append([currentlabel -1, x,y,z,val,float(wl),val/noise_at_wl, N_labeled_new, currentlabel])

    new_regs = np.array(new_regs)

    N = len(c)

    print()
    print("{} out of {} manually labeled objects had associated segments in original map. Added others.".format(count_labeled,N))
    print("Original automatic detections have ids starting with 1.")
    print("Original automatic detections that were also flagged manually have ids starting with {}.".format(auto_and_man_label_start))
    print("Detections that were only flagged manually have ids starting with {}.".format(only_man_label_start))
    print("Max Segment ID is: ", np.max(m))
    
    return m



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sigma_grow_threshold", type=float,  default=1.7,
                        help="Region growth threshold.")
    parser.add_argument('-i', '--infile', type=str,                    
                        help='Input cube.')
    parser.add_argument('-m', '--mapcube', type=str,                    
                        help='Input map cube.')
    parser.add_argument('-o', '--outfile', type=str,
                        help='Output map.')
    parser.add_argument('-c', '--catalog', type=str,
                        help='Manual catalog file.')
    parser.add_argument('-n', '--noisemodel', type=str,
                        help='File containing the polynomial parameters for the noise model.')

    args = parser.parse_args(sys.argv[1:])

    catfile  = args.catalog
    mapcube  = args.mapcube
    datacube = args.infile
    noisemodel = args.noisemodel
    sigma_grow_threshold  = args.sigma_grow_threshold
    newmapcube = args.outfile

    #load noise model
    p = np.loadtxt(noisemodel)

    # read manual catalog
    with open(catfile, 'r') as f:
        ll = f.readlines()
    c = []
    for l in ll:
        if l.strip().startswith("#"):
            continue
        if l.strip() == "":
            continue
        x,y,z = [int(np.round(float(f))) for f in l.split()[:3]]
        c.append([x,y,z])

    # open data and map cubes
    mhdu = fits.open(mapcube)
    dhdu = fits.open(datacube)
    m = mhdu[0].data
    d = dhdu[0].data
    wcs = WCS(dhdu[0].header)

    N_labeled_original = np.sum(m > 0.)

    # add segments for manually labeled regions.
    m = apply_manual_catalog(c, m, d, wcs, p, sigma_grow_threshold)

    # save new map
    mhdu[0].data = m
    mhdu.writeto(newmapcube, overwrite=True)
    print("Wrote {}.".format(newmapcube))

   
if __name__=='__main__':
    main()

