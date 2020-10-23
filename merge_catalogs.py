import sys

import glob
import argparse

from astropy.table import vstack
from astropy.table import vstack, Column
from astropy.io import ascii



parser = argparse.ArgumentParser()
parser.add_argument("-f", "--field", type=str, default='COSMOSA',
                    help="Field to merge catalogs for.")


args = parser.parse_args(sys.argv[1:])




ff = glob.glob("data/msf2outcube_{field}_???.cat".format(field=args.field))
#print(ff)
tt = []
for f in ff:
    ifu = f[-7:-4]
    print(ifu)
    t = ascii.read(f, format="ecsv")
    t.add_column(Column([ifu]*len(t),name='ifu'),0)
    tt.append(t)
 
tall = vstack(tt)
tall.write("data/msf2outcube_{field}_allifu.cat".format(field=args.field),format="ascii.ecsv") 