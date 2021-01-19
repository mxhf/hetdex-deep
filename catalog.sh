#!/bin/bash
z=`xpaget single_cube cube slice wcs`
xy=`xpaget single_cube crosshair`
echo $xy $z $1 >> $2
