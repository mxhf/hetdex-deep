#!/bin/bash
ds9 -title single_cube $1 $1 &

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
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits -$vmax 0
xpaset -p single_cube cmap invert

xpaset -p single_cube tile
xpaset -p single_cube tile grid mode manual
xpaset -p single_cube tile grid layout 2 1

xpaset -p single_cube zoom to fit

xpaset -p single_cube crosshair

xpaset -p single_cube frame 2
xpaset -p single_cube catalog sdss9
xpaset -p single_cube frame 1
xpaset -p single_cube catalog symbol color white
xpaset -p single_cube lock crosshair wcs

#frame 7x4
