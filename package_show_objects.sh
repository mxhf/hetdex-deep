#!/bin/bash

mkdir -p pkg_show_objects
cd pkg_show_objects
cp ../pdz_cosmos2015_v1.3.fits.gz .
mkdir -p data
cd data
ln -s ../../data/mmsf2outcube_COSMOS?_allifu.cat .
ln -s ../../data/sf2outcube_COSMOS?_???.fits.gz .
ln -s ../../data/outcube_COSMOS?_???.fits.gz .
ln -s ../../data/mmap_COSMOS?_???.fits.gz .
ln -s ../../data/msf2outcube_COSMOS?_allifu.cat .
cd ..
mkdir -p aux/thumbnails/COSMOSC 
cd aux/thumbnails/COSMOSC
ln -s ../../../../aux/thumbnails/COSMOSC/ifu???_*.fits .
cd ../../../

mkdir -p aux/thumbnails/COSMOSD
cd aux/thumbnails/COSMOSD
ln -s ../../../../aux/thumbnails/COSMOSD/ifu???_*.fits .
cd ../../../


mkdir -p src
cp ../src/show_objects.ipynb src/.
cp ../src/spectrum.py src/.
cd ..
cd ..
tar -cvf --dereference pkg_show_objects.tar pkg_show_objects
