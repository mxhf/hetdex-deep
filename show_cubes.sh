#!/bin/bash
ds9 data/sf2outcube_median_???.fits.gz &

echo "Hit ENTER when done loading..."
read foo

vmax=0.6
xpaset -p ds9 frame lock image
xpaset -p ds9 cube lock wcs
xpaset -p ds9 frame lock colorbar
xpaset -p ds9 scale lock
xpaset -p ds9 scale lock limits

xpaset -p ds9 cmap staircase
xpaset -p ds9 scale limits 0 $vmax

xpaset -p ds9 tile
xpaset -p ds9 tile grid mode manual
xpaset -p ds9 tile grid layout 7 4

xpaset -p ds9 zoom to fit
xpaset -p ds9 colorbar orientation vertical
#frame 7x4

function save_regs(){
xpaset -p ds9 frame $1; xpaset -p ds9 regions save frame$1.reg
}
save_regs 1
