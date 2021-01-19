#!/bin/bash
wl=`xpaget single_cube cube wcs`
python ../src/auto_limits.py $noisemodel $wl
