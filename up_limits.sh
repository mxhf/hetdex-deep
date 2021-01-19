#!/bin/bash


vmax=`xpaget single_cube scale limits | awk '{print $2+.1}'`
xpaset -p single_cube frame 1
xpaset -p single_cube scale limits 0. $vmax

xpaset -p single_cube frame 2
xpaset -p single_cube scale limits -$vmax 0.

xpaset -p single_cube frame 3
xpaset -p single_cube scale limits -$vmax 0.

xpaset -p single_cube frame 1

