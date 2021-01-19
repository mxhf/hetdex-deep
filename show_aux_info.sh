
#!/bin/bash


fncube=`xpaget single_cube file | sed 's/\[/ /g' | awk '{print $1}'`
field=`xpaget single_cube file | sed 's/\[/ /g' | sed 's/_/ /g' | awk '{print $2}'`
ifu=`xpaget single_cube file | sed 's/\[/ /g' | sed 's/_/ /g' | sed 's/\./ /g' | awk '{print $3}'`
z=`xpaget single_cube cube IMAGE`
xy=`xpaget single_cube crosshair`
vmax=`xpaget single_cube scale limits | awk '{print $2}'`
levels=`awk -v vmax=$vmax 'BEGIN{s=vmax/0.6; print .2*s, .3*s, .4*s, .5*s, .6*s, .8*s, 1.0*s}'`

echo $fncube $field $ifu $xy $z $1

fnthumbB=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.B.original_psf.v2.fits
fnthumbV=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.V.original_psf.v2.fits
fnthumbgp=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.gp.original_psf.v2.fits
fnthumbrp=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.rp.original_psf.v2.fits
fnthumbip=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.ip.original_psf.v2.fits
fnthumbzp=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.zp.original_psf.v2.fits
fnthumbKs=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.Ks.original_psf.v5.fits
#fnthumbK=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.Ks.original_psf.v5.fits
fnthumbK=../aux/thumbnails/${field}/ifu${ifu}_COSMOS.K.UV_original_psf.v1.fits

instance=aux_info
a=`xpaaccess $instance `
b=`xpaget $instance about 2> /dev/null | grep SAOImage | wc | awk '{print $1}'`
c=`xpaget $instance file | grep fits | wc | awk '{print $1}'`

echo $a $b $c
#if [[ ( $a != "yes" || $b == 0 || $c == 0 ]]; then
#if [[ ( $a != "yes" ]]; then
#	echo "Starting ds9."
#fi
ds9 -geometry 1440x700  -threads 6 -title aux_info $fncube $fnthumbB $fnthumbV $fnthumbgp $fnthumbrp $fnthumbip $fnthumbzp $fnthumbKs $fnthumbK&




### wait for ds9 to load
instance=aux_info
a="no"
b=0
n=20
dt=0.1
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

xpaset -p aux_info tile grid mode manual
xpaset -p aux_info tile grid layout 5 2

xpaset -p aux_info frame lock wcs
xpaset -p aux_info cube lock wcs
#xpaset -p aux_info frame lock colorbar
#xpaset -p aux_info scale lock
#xpaset -p aux_info scale lock limits


xpaset -p aux_info frame 1
xpaset -p aux_info pan to $xy image
xpaset -p aux_info cube $z

xpaset -p aux_info frame 1
xpaset -p aux_info cmap staircase
xpaset -p aux_info scale limits 0 $vmax

xpaset -p aux_info crosshair

xpaset -p aux_info frame 1
xpaset -p aux_info catalog sdss9
xpaset -p aux_info frame 1
xpaset -p aux_info catalog symbol color white
xpaset -p aux_info lock crosshair wcs
xpaset -p aux_info catalog close

xpaset -p aux_info catalog load ../civano+2016_cls.xml
xpaset -p aux_info catalog symbol color green
xpaset -p aux_info catalog close

function set_color(){
	xpaset -p aux_info frame $1
	xpaset -p aux_info cmap heat
	xpaset -p aux_info scale mode 98.
}
set_color 2
set_color 3
set_color 4
set_color 5
set_color 6
set_color 7
set_color 8
set_color 9


#xpaset -p aux_info frame 1
#xpaset -p aux_info pan to $xy image
#xpaset -p aux_info cube $z

xpaset -p aux_info frame 1
xpaset -p aux_info contour method smooth
xpaset -p aux_info contour generate
xpaset -p aux_info contour 
xpaset -p aux_info contour levels "{${levels}}"
#xpaset -p aux_info contour levels "{.2 .3 .4 .5 .6 .8 1.0}"
xpaset -p aux_info contour color grey
xpaset -p aux_info contour copy

function paste_contours(){
	xpaset -p aux_info frame $1
	xpaset -p aux_info contour paste wcs white 1 no 
}
paste_contours 2
paste_contours 3
paste_contours 4
paste_contours 5
paste_contours 6
paste_contours 7
paste_contours 8
paste_contours 9
xpaset -p aux_info contours close

file=sf2outcube_${field}_${ifu}_photz.reg
if [ -f "$file" ]
then
	echo "$file found."
else
	echo "$file not found."
	python ../src/zphotds9.py ${fncube}
fi

xpaset -p aux_info frame 2
xpaset -p aux_info regions load ${file}  color grey
xpaset -p aux_info frame 1

