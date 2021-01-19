#!/bin/bash


xpaset -p single_cube frame 1
xpaset -p single_cube scale limits 0. $1

xpaset -p single_cube frame 2
xpaset -p single_cube scale limits -$1 0.

xpaset -p single_cube frame 3
xpaset -p single_cube scale limits -$1 0.

xpaset -p single_cube frame 1

