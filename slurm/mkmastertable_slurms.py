#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse
from glob import glob
from astropy.io import ascii

cmd="python ~/hetdex/hetdex-deep/create_master_table.py {shotid}"

slurmtempl="run_mastertable.slurm"


parser = argparse.ArgumentParser()
parser.add_argument('-N', '--njobs', default=48, type=int)
parser.add_argument('-s', '--slurms_dir',  default="slurms_mastertable", type=str)

args = parser.parse_args()

slurms_dir = args.slurms_dir

if not os.path.exists(args.slurms_dir):
    os.mkdir(args.slurms_dir)

worklist = ascii.read("/work/04287/mxhf/maverick/hetdex/shotlist_hdr1.txt")

N = len(worklist)
ii = np.arange(N)
iii = np.array_split(ii,int(args.njobs))

print("Distributing {} tasks over {} jobs.".format(N,len(iii)))
for j,ii in enumerate(iii):
    runfile="run_mastertable_{:04d}.run".format(j)
    with open(slurmtempl, 'r') as f:
        s = f.read()
    s = s.replace("@@@TASKS@@@", str(len(ii)) )
    s = s.replace("@@@CONTROL_FILE@@@", os.path.join(slurms_dir,runfile))
    with open(os.path.join(slurms_dir,slurmtempl.replace(".slurm", "_{:04d}.slurm".format(j))),'w') as f:
        f.write(s)
    path = os.path.join(slurms_dir,runfile)
    with open(path, 'w') as f:
        for i in ii:
            shotid = worklist["datevobs"][i]
            f.write(cmd.format(shotid = shotid) + "\n")
        print("Writing {}".format( path ) )

