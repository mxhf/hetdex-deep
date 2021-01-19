#!/bin/bash

#ds9 -threads 6 -geometry 1440x700 -title ds9

#echo "Hit ENTER when done loading..."
#read foo

vmax=0.6

xpaset -p ds9 cube lock wcs
xpaset -p ds9 frame lock wcs

xpaset -p ds9 zoom to fit 

xpaset -p ds9 cmap staircase

xpaset -p ds9 colorbar lock
xpaset -p ds9 colorbar staircase

xpaset -p ds9 scale lock limits

xpaset -p ds9 frame 1
xpaset -p ds9 scale limits 0 $vmax

