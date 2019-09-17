#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse


cmd="""python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a LL LU -b LU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a LL LU -b LL --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a RL RU -b RU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa={sky_kappa} --ncomp={ncomp} -i {ifu} -a RL RU -b RU --shotlist_pca {shotlist_pca} --shotlist_skyrecon {shotlist_skyrecon} --dir_rebin {dir_rebin}"""

slurmtempl="run_pca.slurm"


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifulistfile', default="ifulist_COSMOSA.txt", type=str)
parser.add_argument('-r', '--dir_rebin', default="/scratch/04287/mxhf/rebin2", type=str)
parser.add_argument('-p', '--shotlist_pca', default="shotlist_PCA_COSMOSA.txt", type=str)
parser.add_argument('-t', '--shotlist_skyrecon', default="shotlist_COSMOSA.txt", type=str)
parser.add_argument('-k', '--sky_kappa', default=5., type=float)
parser.add_argument('-n', '--ncomp', default=20, type=int)
parser.add_argument('-s', '--slurms_dir', default='slurms', type=str)
parser.add_argument('-N', '--njobs', default=48, type=int)

args = parser.parse_args()

ifulistfile = args.ifulistfile
dir_rebin = args.dir_rebin
shotlist_pca = args.shotlist_pca
shotlist_skyrecon = args.shotlist_skyrecon
sky_kappa = args.sky_kappa
ncomp = args.ncomp
slurms_dir = args.slurms_dir

if not os.path.exists(args.slurms_dir):
    os.mkdir(args.slurms_dir)

ifus = []
def loadlist(listfile):
    li = []
    with open(listfile,'r') as f:
        ll = f.readlines()
    for l in ll:
        li.append(l.strip())
    return li

ifus = loadlist(ifulistfile)

worklist = []
for i in ifus:
    worklist.append(i)

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
            ifu = worklist[i]
            f.write(cmd.format(ifu = ifu, dir_rebin = dir_rebin, shotlist_pca = shotlist_pca, shotlist_skyrecon = shotlist_skyrecon, sky_kappa = sky_kappa, ncomp = ncomp) + "\n")
        print("Writing {}".format( path ) )

