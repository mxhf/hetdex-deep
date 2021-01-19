#!/bin/bash
ds9 -title single_cube data/sf2outcube_median_$1.fits.gz aux/bulk/subaru_mosaics/ifu$1_COSMOS.ip.original_psf.v2.fits data/sf2outcube_median_$1.fits.gz map_$1.fits.gz &

echo "Hit ENTER when done loading..."
read foo

vmax=0.6
xpaset -p single_cube frame lock wcs
xpaset -p single_cube cube lock wcs
#xpaset -p single_cube frame lock colorbar
#xpaset -p single_cube scale lock
#xpaset -p single_cube scale lock limits

xpaset -p single_cube frame 1
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits 0 $vmax

xpaset -p single_cube frame 2
xpaset -p single_cube cmap rainbow
xpaset -p single_cube scale limits 0 80

xpaset -p single_cube frame 3
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits -$vmax 0
xpaset -p single_cube cmap invert

xpaset -p single_cube frame 4
xpaset -p single_cube cmap b
xpaset -p single_cube scale asinh
xpaset -p single_cube scale mode 99.

xpaset -p single_cube tile
xpaset -p single_cube tile grid mode manual
xpaset -p single_cube tile grid layout 2 2

xpaset -p single_cube zoom to fit
#frame 7x4
