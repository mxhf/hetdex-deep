#!/bin/bash
#ds9 -title single_cube $1 $1 $2 $3&

sigma=${3:-2}
echo sigma $sigma

vmax=0.6

export catalog=$1_$2.cat
export noisemodel=sf${sigma}outcube_$1_$2.detect_noise_model

#export catalog=$1_$2.txt
#ds9 -threads 6 -geometry 1440x700 -title single_cube outcube_$1_$2.fits.gz outcube_$1_$2.fits.gz sf2ncoutcube_$1_$2.fits.gz map_$1_$2.fits.gz &
#ds9 -threads 6 -geometry 1440x700 -title single_cube sf2outcube_$1_$2.fits.gz sf2outcube_$1_$2.fits.gz sf2ncoutcube_$1_$2.fits.gz map_$1_$2.fits.gz &
ds9 -threads 6 -geometry 1440x700 -title single_cube sf${sigma}outcube_$1_$2.fits.gz sf${sigma}outcube_$1_$2.fits.gz sf${sigma}ncoutcube_$1_$2.fits.gz map_$1_$2.fits.gz outcube_$1_$2.fits.gz&
### wait for ds9 to load
instance=single_cube
a="no"
b=0
n=20
dt=0.25
while [[ ( $a != "yes" || $b == 0 || $c == 0) && $i != $n ]]
do
        echo "ds9 not acessible yet. Tying `awk -v n=$n -v i=$i 'BEGIN{print n-i}'` more times."
        i=$[$i+1]
        sleep $dt
        a=`xpaaccess $instance `
        b=`xpaget $instance about 2> /dev/null | grep SAOImage | wc | awk '{print $1}'`
        c=`xpaget $instance file | grep fits | wc | awk '{print $1}'`
done
if [[ $i == $n ]];
then
        echo "ds9 after $n tries still not acessible, giving up."
        exit 1
fi
echo "Done loading"
### end wait for ds9 to load
#echo "Hit ENTER when done loading..."
#read foo


echo "Writing to catalog " $catalog

xpaset -p single_cube tile column
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

xpaset -p single_cube frame 3
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits -$vmax 0
xpaset -p single_cube cmap invert

xpaset -p single_cube frame 4
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits -$vmax 0
xpaset -p single_cube cmap invert

xpaset -p single_cube frame 5
xpaset -p single_cube cmap staircase
xpaset -p single_cube scale limits 0 `python -c "print($vmax*2)"`

xpaset -p single_cube tile
xpaset -p single_cube tile grid mode manual
xpaset -p single_cube tile grid layout 4 1

xpaset -p single_cube zoom to fit

xpaset -p single_cube crosshair

xpaset -p single_cube frame 1
xpaset -p single_cube catalog sdss9
xpaset -p single_cube frame 1
xpaset -p single_cube catalog symbol color white
xpaset -p single_cube catalog close

xpaset -p single_cube crosshair 70 70
xpaset -p single_cube lock crosshair wcs


xpaset -p single_cube catalog load ../civano+2016_cls.xml
xpaset -p single_cube catalog symbol color green
xpaset -p single_cube catalog close

#frame 7x4

#xpaset -p single_cube cube play
xpaset -p single_cube cube interval .125

tail -f $catalog

xpaset -p single_cube analysis load ../show_cube3.ans
xpaset -p single_cube analysis load ../show_cube3.ans
