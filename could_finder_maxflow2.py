
# coding: utf-8

# In[2]:

COMMANDLINE = True
TESTING = False


# In[3]:

if not COMMANDLINE:
    # go wide screen
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:100% !important; }</style>"))


# In[4]:

# This version implements a 
# spectrally variing detection threshold.
import time
import numpy as np
import spectrum
from collections import OrderedDict
import argparse
from astropy.io import fits
from astropy import wcs
import sys
from astropy.io import ascii
import os

import numpy as np
#from imageio import imread
if not COMMANDLINE:
    get_ipython().magic('matplotlib inline')
else:
    import matplotlib
    matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import maxflow


# ## fast graph construction, no color or distance based weights

# In[42]:

if TESTING:
    fcube = "sf2outcube_COSMOSA_022_pca.fits.gz"
    fnoisemodel = "sf2outcube_COSMOSA_022.detect_noise_model"
    workdir = "/Users/mxhf/ownCloudRZG/work/MPE/hetdex/src/deep_and_repeat/COSMOS"
    datadir = os.path.join(workdir, "data")

    w = 10.
    wz = .05
    dscale = 2.5

    # read data
    s = spectrum.readSpectrum(os.path.join(datadir, fcube))
    ww = s.grid()

    # read noise model
    noise_model = np.loadtxt(os.path.join(datadir, fnoisemodel) )

    # scale date by noise
    c = s.data

    c = dscale * np.multiply(c.T, 1./np.polyval(noise_model, ww )).T

    c = c[118-5:118+6] # for faster testing




    # Create a graph with integer capacities.
    # best guess for number of nodes and edges
    nnodes = len( c.flatten() )
    nedges = len( c.flatten() )*26 + 2*len( c.flatten() )
    g = maxflow.Graph[float](nnodes, nedges )

    nodes = g.add_grid_nodes(c.shape)

    # connection weights based on Euclidian distances
    structure = np.array([ w*wz* np.array([[.58, .7, .58],
                            [.7, 1, .7],
                            [.58, .7, .58]]),
                           w*np.array([[.7, 1, .7],
                            [1, 0, 1],
                            [.7, 1, .7]]),
                           w*wz*np.array([[.58, .7, .58],
                            [.7, 1, .7],
                            [.58, .7, .58]]) ])

    g.add_grid_edges(nodes, structure=structure)
    g.add_grid_tedges(nodes, c, 1.-c)


    # Find the maximum flow.
    g.maxflow()

    # Get the segments.
    sgm = g.get_grid_segments(nodes)

    # The labels should be 1 where sgm is False and 0 otherwise.
    img2 = np.int_(np.logical_not(sgm))
    # Show the result.
    cout = img2.reshape(c.shape)

if TESTING:
    #subcube = c[118-2:118+3]
    #subcout = cout[118-2:118+3]
    subcube = c
    subcout = cout
    N = subcube.shape[0]
    if N < 15:
        f = plt.figure(figsize=[15,10])
        for i in range(N):
            plt.subplot(1,N,i+1)
            plt.imshow(subcube[i]/dscale, origin='bottom', interpolation='none',vmin=0., vmax=2.)

        f = plt.figure(figsize=[15,10])
        for i in range(N):
            plt.subplot(1,N,i+1)
            plt.imshow(subcout[i], origin='bottom', interpolation='none')
    plt.show()


# # production code

# In[46]:

def maxflow_detect(ww, c, noise_model, dscale, w, wz):
    #print("maxflow_detect: ", ww, c, noise_model, dscale, w, wz)
    """Maxflow/graphcut based detection algorithm."""
    #print("Running maxflow/graphcut")
    #1/0
    # scale date by noise
    sc = dscale * np.multiply(c.T, 1./np.polyval(noise_model, ww )).T

    # Create a graph with integer capacities.
    # best guess for number of nodes and edges
    nnodes = len( sc.flatten() )
    nedges = len( sc.flatten() )*26 + 2*len( sc.flatten() )
    g = maxflow.Graph[float](nnodes, nedges )

    nodes = g.add_grid_nodes(sc.shape)

    # connection weights based on Euclidian distances
    structure = np.array([ w*wz* np.array([[.58, .7, .58],
                            [.7, 1, .7],
                            [.58, .7, .58]]),
                           w*np.array([[.7, 1, .7],
                            [1, 0, 1],
                            [.7, 1, .7]]),
                           w*wz*np.array([[.58, .7, .58],
                            [.7, 1, .7],
                            [.58, .7, .58]]) ])

    g.add_grid_edges(nodes, structure=structure)
    g.add_grid_tedges(nodes, sc, 1.-sc)


    # Find the maximum flow.
    g.maxflow()

    # Get the segments.
    sgm = g.get_grid_segments(nodes)

    # The labels should be 1 where sgm is False and 0 otherwise.
    img2 = np.int_(np.logical_not(sgm))
    
    # Show the result.
    detectseeds = img2.reshape(c.shape)
    return detectseeds == 1


# In[49]:

# This version implements a 
# spectrally variing detection threshold.
import time
import numpy as np
import spectrum
from collections import OrderedDict
import argparse
from astropy.io import fits
from astropy import wcs
import sys
from astropy.io import ascii
import os
import maxflow


def pp(s):
    print(s)
    return s + "\n"


def grow_segment(c, outmap, x, y, z, detectseeds, threshold, label):
    # Non object oriented version
    def get_children(zyx, shape):
        maxz, maxy, maxx = shape
        z, y, x = zyx
        children = []
        ddx = [-1,0,1]
        ddy = [-1,0,1]
        ddz = [-1,0,1]
        #ddz = [0]
        for dx in ddx:
            for dy in ddy:
                for dz in ddz:
                    if not dx == dy == dz == 0:
                        newx, newy, newz = x+dx, y+dy, z+dz
                        if newx >= 0 and newx < maxx                            and newy >= 0 and newy < maxy                            and newz >= 0 and newz < maxz:
                               children.append( (newz, newy, newx) ) 
        return children

    pixel_list = []
    stack = [(z,y,x)]
    outmap[(z,y,x)] = label

    maxiter = 1e6
    iter = 0
    while len(stack) > 0: 
        pm = stack.pop()
        pp = get_children(pm, c.shape)
        for p in pp:
            if ((c[p] >= threshold) or (detectseeds[p])) and  outmap[p] == 0: # pixel that are labeled 0 have not been visited yet
                stack.append(p)
                outmap[p] = label
            elif ((c[p] >= threshold) or (detectseeds[p])) and  outmap[p] == 0:
                outmap[p] = -1 # label pixel with -1 if they have beed visited already
        iter += 1
        if iter > maxiter:
            break

    return outmap#pixel_list            



def build_map(detectseeds, ww, c, sigma_detect_threshold, sigma_grow_threshold, noise_model):
    # run cloud growth algorithm on all pixels seeded pixels
    # try to be intelligent, only loop over pixels that exceed detection threshold


    ii = detectseeds
    
    outmap = np.zeros_like(c, dtype=int)

    zz,yy,xx = [np.arange(s, dtype=int) for s in c.shape]
    YY,ZZ,XX = np.meshgrid(yy, zz, xx)

    label = 1
    for i, (x,y,z) in enumerate( zip(XX[ii], YY[ii], ZZ[ii] ) )  :  
        if outmap[z,y,x] == 0:
            pp("### {} z={} ###".format(i,z))
            print(x,y,z, outmap[z,y,x])
            #outmap = build_map2(c, outmap, x,y,z, threshold = grow_threshold, label=label)

            summary = ""
            summary += pp("Building map starting with pixel {} out of {} that exceeds threshold...".format(i,N))
            start_time = time.time()
            #print("{} labeled pixels on map".format( np.nansum( (outmap > 0).flatten())) )
            
            # Here we set to which threshold a region shoudl be grown.
            # note that this is not stricly correct, the threshold should very also inside
            # grow_segment as function of wavelength. Here we will just rely on this
            # being a slowly variing function.
            threshold = np.polyval(noise_model, ww[z]) * sigma_grow_threshold
            print("grow_threshold = ", threshold)
            outmap = grow_segment(c, outmap, x,y,z, detectseeds, threshold, label=label)
            time_to_build = time.time() - start_time
            summary += pp("Time to build map: {:.4e} s".format(time_to_build))

            print("{} labeled pixels on map".format( np.nansum( (outmap > 0).flatten())) )
            print("{} untouched pixels on map".format( np.nansum( (outmap == 0).flatten())) )

        label += 1

        if i > 1e6:
            break
    #print (cnt)
    return outmap


def filter_minsize(outmap, minsize):
    # filter regions with too small volume 
    rr = np.sort( np.unique( outmap.flatten() ) )
    for r in rr:
        if not r > 0:
            continue
        N = np.sum( outmap == r )
        #print("{} : N = {}".format(r, N))
        if N < minsize:
            outmap[outmap == r] = -1  
    rr = np.sort( np.unique( outmap.flatten() ) )
    
    # relabel
    for i,r in enumerate(rr[rr>0]):
        outmap[outmap == r] = i + 1

    rr = np.sort( np.unique( outmap.flatten() ) )
    print("{} regions survive size cut".format( len(rr[rr>0]) ))

    return outmap



def save_map(outmap, fmapout):
    w = wcs.WCS(s.hdu)
    # save map
    f = np.zeros_like(c)
    f[outmap > 0] = c[outmap > 0]
    wcs_header =  w.to_header()


    h = fits.PrimaryHDU(data=outmap, header=s.hdu.header)
    for k in wcs_header:
        h.header[k] = wcs_header[k]
    hdu = fits.HDUList(h)

    # save map filtered data
    f = np.zeros_like(c)
    f[outmap > 0] = c[outmap > 0]
    h = fits.ImageHDU(data=f, header=s.hdu.header, name = "filtered_data")
    for k in wcs_header:
        h.header[k] = wcs_header[k]

    hdu.append(h)

    # save shells 
    f = np.zeros_like(c)
    f[outmap == -1] = c[outmap == -1]

    h = fits.ImageHDU(data=f, header=s.hdu.header, name = "shells")
    for k in wcs_header:
        h.header[k] = wcs_header[k]
    hdu.append(h)

    hdu.writeto(fmapout, overwrite=True)
    print("Wrote {}.".format(fmapout))


def threshold_detect(ww, c, noise_model, sigma_detect_threshold):
    """Here we simply flag all pixel that have a value of sigma_detect_threshold x noise.""" 
    detectseeds = (c.swapaxes(2,0) > np.polyval(noise_model, ww) * sigma_detect_threshold ).swapaxes(2,0)
    return detectseeds

if COMMANDLINE:

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--minsize", type=int, default=3,
                        help="Minimum region size.")
    parser.add_argument("-g", "--sigma_grow_threshold", type=float,  default=3.,
                        help="Region growth threshold.")
    parser.add_argument("-d","--sigma_detect_threshold", type=float, default=5.,
                        help="Detection threshold.")
    parser.add_argument('-i', '--infile', type=str,                    
                        help='Input cube.')
    parser.add_argument('-o', '--outfile', type=str, default='',
                        help='Output map.')
    parser.add_argument('-b', '--bad_regions', type=str, default='',
                        help='Output map.')

    parser.add_argument('-n', '--noisemodel', type=str, default='',
                        help='File containing the polynomial parameters for the noise model.')
    
    parser.add_argument('-w', type=float, default=10.,
                        help='Maxflow spatial weight.')
    parser.add_argument('-z',  type=float, default=0.05,
                        help='Maxflow dispersion weight.')
    parser.add_argument('-s', '--dscale', type=float, default=2.5,
                        help='Maxflow datascale.')

    args = parser.parse_args(sys.argv[1:])


    fcube = args.infile
    fmapout = args.outfile
    sigma_detect_threshold = args.sigma_detect_threshold
    sigma_grow_threshold = args.sigma_grow_threshold
    minsize = args.minsize
    noise_model = np.loadtxt(args.noisemodel)
    bad_regions = args.bad_regions
    w = args.w
    wz = args.z
    dscale = args.dscale
    
else:
    
    workdir = "/Users/mxhf/ownCloudRZG/work/MPE/hetdex/src/deep_and_repeat/COSMOS/data"
    fcube = workdir + "/sf2outcube_COSMOSA_022_pca.fits.gz"
    fnoisemodel = workdir + "/sf2outcube_COSMOSA_022.detect_noise_model"
    fmapout = workdir + "/map_COSMOSA_022_pca.fits.gz"

    
    bad_regions = ""
    sigma_detect_threshold = 5.
    sigma_grow_threshold = 3.
    minsize = 3.
    print("Loading noisemodel: ", fnoisemodel )
    noise_model = np.loadtxt(fnoisemodel)


    w = 10.
    wz = .05
    dscale = 2.5


# In[51]:

badwlregions = []

if bad_regions != '':
    if not os.path.exists(args.bad_regions):
        print("WARNING: {} not found proceeding without.".format(args.bad_regions))
    else:
        t = ascii.read(args.bad_regions)
        for r in t:
            badwlregions.append([r["start"],r["end"]])

print("Reading input datacube {}.".format(fcube))   
s = spectrum.readSpectrum(fcube)
ww = s.grid()
c = s.data
for r in badwlregions:  
    c[r[0]:r[1]] = 0. # take out bad region
    
# The detects seeds array is a boolean array of the same dimensions
# as the input datacube.
# It flags pixels as potential objects.
detectseeds1 = threshold_detect(ww, c, noise_model, sigma_detect_threshold)


    
# now use maxflow to detect low surface brightness structures.  
detectseeds2 = maxflow_detect(ww, c, noise_model, dscale, w, wz)


# In[52]:

detectseeds = detectseeds1 + detectseeds2
#detectseeds = detectseeds2

N =  (np.sum(detectseeds))
print("{} pixel are flagged as detection seeds".format(N))



outmap = build_map(detectseeds, ww, c, sigma_detect_threshold, sigma_grow_threshold, noise_model)
outmap = filter_minsize(outmap, minsize)


save_map(outmap, fmapout)


# In[ ]:

# python src/could_finder_maxflow2.py --sigma_grow_threshold 3. --sigma_detect_threshold 5. --infile data/sf2outcube_COSMOSA_022_pca.fits.gz -n data/sf2outcube_COSMOSA_022.detect_noise_model -o data/map_COSMOSA_022_pca.fits.gz


# In[ ]:



