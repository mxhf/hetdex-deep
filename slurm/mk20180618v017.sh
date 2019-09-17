#/bin/bash

ifu=024
shot=20180618v017
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa=1.5 --ncomp=25 -i $ifu -a LL LU -b LU --shotlist_pca shotlist_PCA_${shot}.txt --shotlist_skyrecon shotlist_${shot}.txt --dir_rebin /scratch/04287/mxhf/rebin2
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa=1.5 --ncomp=25 -i $ifu -a LL LU -b LL --shotlist_pca shotlist_PCA_${shot}.txt --shotlist_skyrecon shotlist_${shot}.txt --dir_rebin /scratch/04287/mxhf/rebin2
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa=1.5 --ncomp=25 -i $ifu -a RL RU -b RU --shotlist_pca shotlist_PCA_${shot}.txt --shotlist_skyrecon shotlist_${shot}.txt --dir_rebin /scratch/04287/mxhf/rebin2
python /home1/04287/mxhf/hetdex/xpca/pca_sky6.py --sky_kappa=1.5 --ncomp=25 -i $ifu -a RL RU -b RL --shotlist_pca shotlist_PCA_${shot}.txt --shotlist_skyrecon shotlist_${shot}.txt --dir_rebin /scratch/04287/mxhf/rebin2
python ~/hetdex/hetdex_cube/global_sky2.py --dir_rebin=/scratch/04287/mxhf/rebin2  --shot=${shot}
python ~/hetdex/hetdex_cube/cube.py --no_pca --shiftsdir karls_shifts  --dir_rebin=/scratch/04287/mxhf/rebin2 --ifuslot ${ifu} --shotlist shotlist_${shot}.txt -o outcube_${shot}_${ifu}_nopca.fits.gz
python ~/hetdex/hetdex_cube/cube.py --shiftsdir karls_shifts  --dir_rebin=/scratch/04287/mxhf/rebin2 --ifuslot ${ifu} --shotlist shotlist_${shot}.txt -o outcube_${shot}_${ifu}.fits.gz
