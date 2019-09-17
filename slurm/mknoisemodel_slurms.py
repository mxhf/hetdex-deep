#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse
from glob import glob

cmd="python ~/hetdex/hetdex_cube/noise_model2.py --infile {filename} --outfile {outfilename}"

slurmtempl="run_noisemodel.slurm"


parser = argparse.ArgumentParser()
parser.add_argument('-N', '--njobs', default=24, type=int)
parser.add_argument('-s', '--slurms_dir',  default="slurms_noisemodel", type=str)

args = parser.parse_args()

slurms_dir = args.slurms_dir

if not os.path.exists(args.slurms_dir):
    os.mkdir(args.slurms_dir)

worklist = glob("data/outcube*.fits.gz")

N = len(worklist)
ii = np.arange(N)
iii = np.array_split(ii,int(args.njobs))

print("Distributing {} ifus over {} jobs.".format(N,len(iii)))
for j,ii in enumerate(iii):
    runfile="run_filter_{:04d}.run".format(j)
    with open(slurmtempl, 'r') as f:
        s = f.read()
    s = s.replace("@@@TASKS@@@", str(len(ii)*4) )
    s = s.replace("@@@CONTROL_FILE@@@", os.path.join(slurms_dir,runfile))
    with open(os.path.join(slurms_dir,slurmtempl.replace(".slurm", "_{:04d}.slurm".format(j))),'w') as f:
        f.write(s)
    path = os.path.join(slurms_dir,runfile)
    with open(path, 'w') as f:
        for i in ii:
            filename = worklist[i]
            h,t = os.path.split(filename)
            outfilename  = os.path.join(h,"nc" + t)
            f.write(cmd.format(filename = filename, outfilename=outfilename) + "\n")
        print("Writing {}".format( path ) )

