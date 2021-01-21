from astropy.io import fits
import spectrum
from matplotlib import pyplot as plt
import numpy as np
from astropy.stats import biweight_midvariance as bwtmv
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--incube", type=str, 
                    help="Input data cube.")
parser.add_argument("-n", "--noisemodel", type=str, 
                    help="Detection noise model.")
parser.add_argument("-o", "--outfile", type=str,
                    help="Input data cube.")


args = parser.parse_args(sys.argv[1:])

print("Reading {} ...".format(args.incube))
s = spectrum.readSpectrum(args.incube)

print("Reading noisemodel {} ...".format(args.noisemodel))
with open(args.noisemodel) as f:
    ll = f.readlines()
p = [float(l) for l in ll]

wl = s.grid()

nn = np.polyval(p,wl)

s.data = (s.data.T/nn).T

s.hdu.writeto(args.outfile, overwrite=True)

print("Wrote {}.".format(args.outfile))
