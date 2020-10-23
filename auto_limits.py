#!/usr/bin/env python
import numpy as np
import sys
import os

nm = sys.argv[1]
wl = float(sys.argv[2])
scale = 5.

p = np.loadtxt(nm)

noiselevel = np.polyval(p, wl)
vmax = noiselevel * scale
print("Noise level is {noiselevel}. Setting vmax = {vmax}.".format(noiselevel=noiselevel,vmax=vmax))

cmd = """
xpaset -p single_cube frame 1
xpaset -p single_cube scale limits 0 {vmax}

xpaset -p single_cube frame 2
xpaset -p single_cube scale limits -{vmax} 0

xpaset -p single_cube frame 3
xpaset -p single_cube scale limits -{vmax} 0

xpaset -p single_cube frame 4
xpaset -p single_cube scale limits -{vmax} 0

xpaset -p single_cube frame 5
xpaset -p single_cube scale limits 0 {vmax2}

xpaset -p single_cube frame 1
""".format(vmax=vmax,vmax2=vmax*2.)
os.system(cmd)
