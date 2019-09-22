import sys
sys.path.append('/home1/04287/mxhf/HETDEX_API/')

import sys
import os
import os.path
import subprocess
import numpy as np

import tables as tb
from astropy.table import Table

from hetdex_api.survey import Survey

shotid = sys.argv[1]

oufilename = "/work/04287/mxhf/maverick/hetdex/cubes{shotid}.images.txt".format(shotid=shotid)
if os.path.exists(oufilename):
    print("{} exists already.".format(oufilename))
    sys.exit(0)

print("Processing shot {} .".format(shotid))
fileh = tb.open_file('/work/03946/hetdex/hdr1/reduction/data/{shotid}.h5'.format(shotid=shotid))
l = Table( fileh.root.Data.Images[:] )
l2 = Table( [l["contid"], l["ifuid"], l["ifuslot"], l["specid"]] )
l2.write(oufilename, format="ascii")

