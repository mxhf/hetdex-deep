#!/bin/bash
#$ -cwd

python ../src/hetdex_cube/cube.py --shiftsdir ../shifts --ifuslot 022 --shotlist shotlist_COSMOSA.txt -o outcube_COSMOSA_022.fits
