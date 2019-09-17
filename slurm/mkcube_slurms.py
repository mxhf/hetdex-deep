#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse


cmd="python ~/hetdex/hetdex_cube/cube.py --no_pca --shiftsdir karls_shifts  --dir_rebin={dir_rebin} --ifuslot {ifu} --shotlist shotlist_{field}.txt -o outcube_{field}_{ifu}_nopca.fits.gz"
slurmtempl="run_cube.slurm"


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--worklist',  type=str)
parser.add_argument('-r', '--dir_rebin', default="/scratch/04287/mxhf/rebin2", type=str)
parser.add_argument('-N', '--njobs', default=48, type=int)
parser.add_argument('-s', '--slurms_dir',  default="slurms_cubes", type=str)

args = parser.parse_args()

worklistfile = args.worklist
dir_rebin = args.dir_rebin
slurms_dir = args.slurms_dir

if not os.path.exists(args.slurms_dir):
    os.mkdir(args.slurms_dir)

def loadlist(listfile):
    li = []
    with open(listfile,'r') as f:
        ll = f.readlines()
    for l in ll:
        li.append(l.strip().split())
    return li

worklist = loadlist(worklistfile)

N = len(worklist)
ii = np.arange(N)
iii = np.array_split(ii,int(args.njobs))

print("Distributing {} ifus over {} jobs.".format(N,len(iii)))
for j,ii in enumerate(iii):
    runfile="run_pca_{:04d}.run".format(j)
    with open(slurmtempl, 'r') as f:
        s = f.read()
    s = s.replace("@@@TASKS@@@", str(len(ii)*4) )
    s = s.replace("@@@CONTROL_FILE@@@", os.path.join(slurms_dir,runfile))
    with open(os.path.join(slurms_dir,slurmtempl.replace(".slurm", "_{:04d}.slurm".format(j))),'w') as f:
        f.write(s)
    path = os.path.join(slurms_dir,runfile)
    with open(path, 'w') as f:
        for i in ii:
            field, ifu = worklist[i]
            f.write(cmd.format(field = field, ifu = ifu, dir_rebin = dir_rebin) + "\n")
        print("Writing {}".format( path ) )

