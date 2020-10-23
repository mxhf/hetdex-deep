USEDS9 = False

if USEDS9:
    import pysao
from astropy.io import ascii, fits
from astropy.table import Table, Column
import os

import ipywidgets as widgets
from IPython.display import display
from ipywidgets import Button, Layout

import numpy as np

from astropy.visualization import (MinMaxInterval, LogStretch,
                           ImageNormalize)

from astropy.visualization import PercentileInterval

from matplotlib import cm

from astropy.wcs import wcs

import spectrum


import matplotlib.pyplot as plt
import numpy as np


from ipywidgets import Layout, widgets
from IPython.display import display,clear_output

vmax = 0.6
imgscale = 1.0
x,y,z = 0.,0.,0.

fncube = ""
fnorigcube = ""
fnmap = ""
fncatalog = ""

def show_next_record():
    global ii, catalog, current_record_idx, info, phz_zz, phz_pdz,c,m,x,y,z
    records = catalog[ii]
    record = records[current_record_idx]
    #print(record["comment2"])
    
    s = ""
    s = "<style>"
    s = "table, th, td {"
    s = "  border: 1px solid black;"
    s = "  border-collapse: collapse;"
    s = "}"
    s = "th, td {"
    s = "  padding: 20px;"
    s = "}"
    s = "</style>"

    s += "<b>Object {} of {}</b><br>\n".format(current_record_idx+1, len(records))
    
    
    ncols = 3
    tablerows = []
    sr = ""
    z_lya = record["wl_com"] /1216. - 1.
    z_oii = record["wl_com"] /3727. - 1.
    
    s += "<b>z_LyA = {:.2f}</b><br>\n".format(z_lya)
    s += "<b>z_OII = {:.2f}</b><br>\n".format(z_oii)
    
    s += "<table style='width:100%'>\n"
    
    for i,n in enumerate(record.colnames):
        sr += "<td><b>{:10s}:</b> ".format(n)
        sr += str(record[n]) + " "
        if record.columns[n].unit != "":
            ustr = str(record.columns[n].unit)
            if ustr == "None":
                ustr = ""
            sr += ustr + "</td>\n"
        else:
            sr += "</td>\n"
        if i % ncols == 0.:
            tablerows.append(sr)
            sr = ""
    s += "</tr>\n"
    s += "</tr>\n<tr>".join(tablerows)       
    s += "</tr>\n"
    s += "</table>\n"
    
    info.value = s
    comment.value = str(record["comment2"])
    
    if current_record_idx == len(records):
        print("Last objects reached.")
        return
    
    
    id = record["id"]
    x = record["manualx"]
    y = record["manualy"]
    z = record["manualz"]
    if np.isnan(x):
        x = record["x_com"]
    if np.isnan(y):
        y = record["y_com"]
    if np.isnan(z):
        z = record["z_com"]
        
    show_object(x,y,z,vmax,gbimages)
    current_record_idx += 1
    
    zlya = record["wl_com"]/1216. - 1.
    zoii = record["wl_com"]/3727. - 1.
    plot_phz(record["ra_com"],record["dec_com"],phz_zz, phz_pdz, zlya, zoii)
    

    plot_spec(c,m,record["id"], record["wl_com"])
    #try:
    #    plot_spec(c,m,record["id"], record["wl_com"])
    #except:
    #    print("Error plotting spectrum.")


def update_catalog(classifier):
    global ii, catalog, current_record_idx, info
    records = catalog[ii]
    
    record = records[current_record_idx-1]
    id = record["id"]
    ifu = record["ifu"]
    jj  = catalog["id"] == id
    jj *= catalog["ifu"] == ifu
    if not np.sum(jj) == 1:
        print("Unable to find object {} for ifu {}.".format(id, ifu))
    catalog["class2"][jj] = classifier
 
def update_comment(comment):
    global ii, catalog, current_record_idx, info
    records = catalog[ii]
    
    record = records[current_record_idx-1]
    id = record["id"]
    ifu = record["ifu"]
    jj  = catalog["id"] == id
    jj *= catalog["ifu"] == ifu
    if not np.sum(jj) == 1:
        print("Unable to find object {} for ifu {}.".format(id, ifu))
    print("@@@ comment = ",comment )
    catalog["comment2"][jj] = comment
 

def on_scale_up(b):
    global vmax
    vmax = float(ds9.xpa.get("scale limits").split()[1])
    ds9.xpa.set("frame 1")
    ds9.xpa.set("scale limits 0. {}".format(vmax*1.05))

def on_scale_down(b):
    global vmax
    vmax = float(ds9.xpa.get("scale limits").split()[1])
    ds9.xpa.set("frame 1")
    ds9.xpa.set("scale limits 0. {}".format(vmax*.95))

    
def on_prevslice_clicked(b):
    ds9.xpa.set("frame 1")
    ds9.xpa.set("cube prev")

def on_nextslice_clicked(b):
    ds9.xpa.set("frame 1")
    ds9.xpa.set("cube next")
    
    
def on_prev_clicked(b):
    global current_record_idx
    current_record_idx -=2
    show_next_record()
    #print("Previous object.")

def on_next_clicked(b):
    global current_record_idx, records
    if current_record_idx == len(records):
        print("No more objects.")
    else:
        show_next_record()
        # always save catalog
        on_save_clicked(None)

    
def on_jump_clicked(b):
    global current_record_idx, records
    idxx = np.arange(len(records))
    tt = object_selection.value.split()
    
    ii = records["id"] == int(tt[0])
    current_record_idx = idxx[ii][0] 
    on_next_clicked(b)
    
def on_halo_clicked(b):
    update_catalog("halo")
    #show_next_record()
    #print("Selected halo.")
    
def on_fil_clicked(b):
    update_catalog("filament")
    #show_next_record()
    #print("Selected halo.")
    
def on_ps_clicked(b):
    update_catalog("pointsource")
    #show_next_record()
    #print("Selected poinsource.")
    
def on_lae_clicked(b):
    update_catalog("LAE")
    
def on_oii_clicked(b):
    update_catalog("OII")
    
def on_gal_clicked(b):
    update_catalog("gal")
    #show_next_record()
    #print("Selected poinsource.")
    
def on_agn_clicked(b):
    update_catalog("agn")
    #show_next_record()
    #print("Selected poinsource.")
    

def on_junk_clicked(b):
    update_catalog("junk")
    #show_next_record()
    #print("Selected junk.")
    
    
def on_save_clicked(b):
    global fncatalog
    #catalog.write(newfncatalog, format="ascii.ecsv", overwrite=True)
    #fncatalog="../data/msf2outcube_{}_allifu.cat".format(field)
    h,t = os.path.split(fncatalog)
    newfncatalog = os.path.join( h , "m" + t)
    catalog.write(newfncatalog, format="ascii.ecsv", overwrite=True)

    
def on_transfer_contours(b):
    transfer_contours()
    
    
def on_clear_contours(b):
    clear_contours_on_images()
    
    
def on_add_comment_clicked(b):
    global comment
    update_comment(comment.value)
    
    
def on_change_ifu(b):
    pass
    
    
if USEDS9:
    def startup(field, ifu, fncube, fnorigcube, fnmap, gbimages, vmax):
        global ds9, LOADSDSS
        ds9 = pysao.ds9()
        ds9.xpa.set("frame 1")
        ds9.xpa.set("fits {fncube}".format(fncube=fncube))

        ds9.xpa.set("frame 2")
        ds9.xpa.set("fits {fnorigcube}".format(fnorigcube=fnorigcube))

        ds9.xpa.set("frame lock wcs")
        ds9.xpa.set("cube lock wcs")

        ds9.xpa.set("tile")
        ds9.xpa.set("tile grid mode manual")
        ds9.xpa.set("tile grid layout 6 2")

        ds9.xpa.set("frame 3")
        ds9.xpa.set("fits {fnmap}".format(fnmap=fnmap))
        for i,f in enumerate(gbimages):
            ds9.xpa.set("frame {}".format(i+4))
            ds9.xpa.set("fits {f}".format(f=f))


        ds9.xpa.set("frame 1")
        ds9.xpa.set("cmap staircase")
        ds9.xpa.set("scale limits 0 {}".format(vmax))
        #ds9.xpa.set("crosshair")

        ds9.xpa.set("frame 2")
        ds9.xpa.set("cmap staircase")
        ds9.xpa.set("scale limits 0 {}".format(vmax*3.))
        #ds9.xpa.set("crosshair")

        ds9.xpa.set("frame 1")
        if LOADSDSS:
            ds9.xpa.set("catalog sdss9")
            ds9.xpa.set("frame 1")
            ds9.xpa.set("catalog symbol color white")
            ds9.xpa.set("catalog close")

        ds9.xpa.set("lock crosshair wcs")
        ds9.xpa.set("catalog load ../civano+2016_cls.xml")
        ds9.xpa.set("catalog symbol color green")
        ds9.xpa.set("catalog close")

        for i,f in enumerate(gbimages):
            ds9.xpa.set("frame {}".format(i+4))
            #ds9.xpa.set("cmap heat")
            ds9.xpa.set("cmap grey")
            ds9.xpa.set("cmap invert yes")
            ds9.xpa.set("scale mode 95.")
            #ds9.xpa.set("scale mode 98.")

        regfile="../data/sf2outcube_{field}_{ifu}_photz.reg".format(field=field,ifu=ifu)
        if not os.path.exists(regfile):
            cmd="python ../src/zphotds9.py {fncube}".format(fncube=fncube)
            #!$cmd
            os.system(cmd)

        ds9.xpa.set("frame 4")
        #print("regions load {}  color grey".format(regfile))
        ds9.xpa.set("regions load {}".format(regfile))
        ds9.xpa.set("regions color grey ")


        ds9.xpa.set("frame 1")

        return ds9

    showing_contours = False


    def clear_contours_on_images():
        global ds9,gbimages, showing_contours
        for i in range(len(gbimages) + 2):
           ds9.xpa.set("frame {}".format(i+1))
           ds9.xpa.set("contour clear ")


    def show_object(x,y,z,vmax,gbimages):
        global showing_contours

        ds9.xpa.set("frame 1")
        ds9.xpa.set("contour clear ")

        if showing_contours:
            clear_contours_on_images()
            showing_contours = False

        ds9.xpa.set("frame 1")
        ds9.xpa.set("pan to {} {} image".format(x,y))
        ds9.xpa.set("crosshair lock wcs")
        ds9.xpa.set("crosshair {} {} image".format(x,y))

        print( "cube {} image".format(int(z)) )
        try:
            # somtimes, realy this is failing, no idea why.
            ds9.xpa.set("cube {} image".format(int(z)))
        except:
            ds9.xpa.set("cube {} image".format(int(z)+1))

        s=vmax/0.6; 
        levels=[ .2*s, .3*s, .4*s, .5*s, .6*s, .8*s, 1.0*s ]
        ds9.xpa.set("frame 1")
        ds9.xpa.set("contour method smooth")
        ds9.xpa.set("contour generate")
        ds9.xpa.set("contour ")
        #print("contour levels " + " ".join([str(l) for l in levels]))
        ds9.xpa.set("contour levels " + " ".join([str(l) for l in levels]))
        ds9.xpa.set("contour color grey")

        ds9.xpa.set("contour close ")

    def transfer_contours():
        global ds9, showing_contours

        ds9.xpa.set("frame 1")
        ds9.xpa.set("contour copy")

        for i,f in enumerate(gbimages):
           ds9.xpa.set("frame {}".format(i+3))
           ds9.xpa.set("contour paste wcs white 1 no ")

        #ds9.xpa.set("contours close")
        ds9.xpa.set("frame 1")
        showing_contours = True


def plot_phz(ra, dec, zz, phz_pdz, zlya, zoii, rmatch = 2.5):
    global phz_canvas
    dd = np.sqrt(((phz_pdz[1].data["RA"] - ra)*np.cos(np.deg2rad(dec)))**2. + (phz_pdz[1].data["DEC"] - dec)**2.)*3600.
    ii = dd < rmatch
    N = np.sum(ii)

    #print("{} objects within {} arcsec.".format( N, rmatch ) )
    
    plt.ioff()
    #ax=plt.gca()
    #ax.clear()
    fig = plt.figure(figsize=[3,3])
    ax = plt.subplot()
        
    if N > 0:
        indices = np.arange(len(phz_pdz[1].data))
        for j in range(N):
            i = indices[ii][j]
            pp = phz_pdz[1].data[i][75:675]
            ax.plot(phz_zz, pp/np.sum(pp), '-')
    else:
        ax.text(.5,.5,"no match",transform=ax.transAxes, ha='center', va='center')
        
    ax.set_xlabel("z")
    ax.set_ylabel("p(z)")
    
    ax.axvline(zlya,c='r')
    ax.axvline(zoii,c='g')
    with phz_canvas:
        clear_output(wait=True)
        display(ax.figure)

        
def plot_spec(c, m, id, wl, win=300.):
    global spec_canvas, current_record_idx

    plt.ioff()
    #ax=plt.gca()
    #ax.clear()
    fig = plt.figure(figsize=[15,3])
    ax = plt.subplot()

    mask = np.sum(  m.data == id, axis=0) > 0
    try:
        spec = np.array( [np.sum(c.data[i][mask]) for i in range(len(c.data))] )
    except:
        spec = np.zeros_like(c.grid())
    
    jj = (c.grid() > (wl - win/2.)) * (c.grid() <= (wl + win/2.))
    vmin,vmax = np.min(spec[jj]),np.max(spec[jj])
    vmin,vmax = -np.max(spec[jj]),np.max(spec[jj])

    
    ax.plot(c.grid(), spec, 'k-')

    ax.set_xlabel("wavelength [A]")
    ax.set_ylabel("flux [arb]")
    
    ax.axvline(wl,c='r')
    ax.set_ylim([vmin,vmax])

    with spec_canvas:
        clear_output(wait=True)
        display(ax.figure)
        
    record = records[current_record_idx]
    # @@@
    p = fit_peak( c.grid(), spec,  record)
    
    z = p[1]/1216. - 1.
    s = "sigma = {:.1f} km/s".format(  p[2]/p[1] * 3e5 ) 
    
    ax.text(.5, .5, s, transform=ax.transAxes, size=20)
    print(s)
    
        
    
    
#ra = 150.22773
#dec =   2.385096
#plot_phz(ra,dec,phz_zz,phz_pdz, zlya=2., zoii=1., rmatch = 2.)

from scipy.optimize import least_squares

import numpy as np


def gauss(mu, sigma, x):
    return 1./(sigma * np.sqrt(2. * np.pi) ) * np.exp( -(x-mu)**2./(2. * sigma**2.))

def peval(p,x):
    A,mu,sigma = p
    return A*gauss(mu, sigma, x)

def resid(p, x, y, yerr=[]):
    model = peval(p,x)
    if yerr==[]:
        return (y - model)
    else:
        return (y - model)/yerr

def fit_gaussians(lineset, ww, csout, wlwin, pp=[]):
    results = []

    for i,wlc in enumerate(lineset):
        if pp != []:
                p0 = pp[i]
        else:
            p0 = [2000.,wlc,10.]
        ii = (ww > wlc-wlwin/2.) * (ww < wlc+wlwin/2.)
        ii *= ~ np.isnan(csout)
        fit = least_squares(resid, p0, args=(ww[ii], csout[ii]))
        p = fit.x

        if False:
            f=plt.figure()

            plt.plot(ww[ii], csout[ii])
            plt.plot(ww[ii], peval(p,ww[ii]))

        results.append([p[0], p[1],p[2]])

    results = np.array(results)


    return results

def fit_peak( ww, csout_unsmoothed,  r, win=300. ):
    print("ww : ", ww)
    ff = csout_unsmoothed
    #yerr = np.std(nsout,axis=0)
    ii = (ww > (r["wl_com"] - win/2.)) * (ww < (r["wl_com"] + win/2.))
    ii *= ~ np.isnan(csout_unsmoothed)

    dw = ww[1]-ww[0]

    A0 = np.sum(ff[ii]) * dw
    mu0    = r["wl_com"]
    sigma0 = r["dwl"]
    p0=[A0,mu0,sigma0]

    #fit = least_squares(resid, p0, args=(ww[ii],ff[ii],yerr[ii]))
    fit = least_squares(resid, p0, args=(ww[ii],ff[ii]))
    p = fit.x
    #print(fit.x)
    ii = (ww > (p[1] - win/2.)) * (ww < (p[1] + win/2.))
    ii *= ~ np.isnan(ff)
    # refit with better centroid
    fit = least_squares(resid, p0, args=(ww[ii],ff[ii]))
    p = fit.x
    #print(fit.x)

    if False:
        f = plt.figure()
        plt.plot(ww[ii],ff[ii])
        plt.errorbar(ww[ii],ff[ii],yerr[ii])
        plt.plot(ww[ii],peval(p,ww[ii]))
        #ax1.fill_between(ww, yerr, -yerr,alpha=0.2, edgecolor='black', facecolor='grey')
    
    return p


def mkgui():
    global info, comment, phz_canvas, spec_canvas, object_selection, figs
    
    button_zoomin = widgets.Button(
        description='zoom in',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_zoomout = widgets.Button(
        description='zoom out',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )


    button_scaleup = widgets.Button(
        description='scale up',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_scaledown = widgets.Button(
        description='scale down',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )


    button_prevslice = widgets.Button(
        description='previous slice',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_nextslice = widgets.Button(
        description='next slice',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_transfer_contours = widgets.Button(
        description='add contours',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_clear_contours = widgets.Button(
        description='clear contours',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_prev = widgets.Button(
        description='Previous',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_next = widgets.Button(
        description='Next',
        disabled=False,
        button_style='success', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_pointsource = widgets.Button(
        description='Pointsource',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_lae = widgets.Button(
        description='LAE',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_oii = widgets.Button(
        description='OII',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_halo = widgets.Button(
        description='Halo',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_fil = widgets.Button(
        description='Filament',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_junk = widgets.Button(
        description='Junk',
        disabled=False,
        button_style='danger', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_gal = widgets.Button(
        description='Galaxy',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_agn = widgets.Button(
        description='AGN',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    button_save = widgets.Button(
        description='Save',
        disabled=False,
        button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )


    info = widgets.HTML(
        value="Hello <b>World</b>",
        placeholder='Some HTML',
        description='Some HTML',
    )


    comment = widgets.Textarea(
        value='',
        placeholder='Type something',
        description='String:',
        disabled=False
    )

    button_add_comment = widgets.Button(
        description='Add comment',
        disabled=False,
        button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='check'
    )

    object_selection = widgets.Dropdown(
        options=['-1'],
        value='-1',
        description='Jump to:',
        disabled=False,
    )
    button_jump = widgets.Button(
        description='Jump',
        disabled=False,
        button_style='success', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Click me',
        icon='arrow'
    )

    button_prevslice.on_click(on_prevslice_clicked)
    button_nextslice.on_click(on_nextslice_clicked)
    button_scaleup.on_click(on_scale_up)
    button_scaledown.on_click(on_scale_down)
    button_transfer_contours.on_click(on_transfer_contours)
    button_clear_contours.on_click(on_clear_contours)
    button_zoomin.on_click(on_zoomin_clicked)
    button_zoomout.on_click(on_zoomout_clicked)

    button_prev.on_click(on_prev_clicked)
    button_next.on_click(on_next_clicked)
    button_pointsource.on_click(on_ps_clicked)
    button_lae.on_click(on_lae_clicked)
    button_oii.on_click(on_oii_clicked)
    button_halo.on_click(on_halo_clicked)
    button_fil.on_click(on_fil_clicked)
    button_gal.on_click(on_gal_clicked)
    button_agn.on_click(on_agn_clicked)
    button_junk.on_click(on_junk_clicked)
    
    button_jump.on_click(on_jump_clicked)

    button_save.on_click(on_save_clicked)
    button_add_comment.on_click(on_add_comment_clicked)

    # buttons, spectra and phz's
    items = [widgets.Label(value="Objects"), button_prev, button_next, \
             button_pointsource, button_lae, button_oii, button_halo, button_fil, button_gal, button_agn, button_junk]
    buttonbox1 = widgets.VBox(items, layout=Layout(width = '150px', positioning="bottom"))

    items = [widgets.Label(value="ds9 control"), button_prevslice, button_nextslice, button_scaleup, \
             button_scaledown, button_zoomin, button_zoomout, button_transfer_contours, button_clear_contours,
             widgets.Label(value="catalog"), button_save]
    buttonbox2 = widgets.VBox(items, layout=Layout(width = '150px', positioning="bottom"))

    
    
    hb5 = widgets.HBox([object_selection, button_jump], layout=Layout(width = '300px', positioning="bottom"))
    
    hb4 = widgets.HBox([buttonbox1, buttonbox2], layout=Layout(width = '300px', positioning="bottom"))
    vb1 = widgets.VBox([hb4, hb5] , layout=Layout(width = '300px', positioning="bottom"))
    
    phz_canvas  = widgets.Output(layout=Layout(height='200px', width = '200px', border='light'))
    spec_canvas = widgets.Output(layout=Layout(height='200px', width = '700px', border='light'))
        
    items = [vb1, info]
    hb1 = widgets.HBox(items, layout=Layout(width='100%', positioning="bottom"))

    items = [comment, button_add_comment]
    hb2 = widgets.HBox(items, layout=Layout(width='100%', positioning="bottom"))
    
    
    hb3 = widgets.HBox([spec_canvas, phz_canvas], layout=Layout(width='100%', positioning="bottom"))

    plt.ioff()
    ax=plt.gca()


    
    if USEDS9:
        return widgets.VBox([hb3, hb1, hb2])
    else:
        # if not using ds9 to display figures, we need additional canvases
        figsize = "200px"

        figs = []
        for i in range(12):
                figs.append( widgets.Output(layout=Layout(height=figsize, width = figsize, border='light')) )

        figures = widgets.VBox([ widgets.HBox(figs[0:4]), widgets.HBox(figs[4:8]), widgets.HBox(figs[8:12])  ])
        return widgets.HBox( [figures, widgets.VBox([hb3, hb1, hb2])] )
    


    object_selection.options = [ "{} {}".format(id, cls2) for id,cls2 in zip( records["id"],records["class2"] ) ]
    

    
def set_input_files(dataroot, field, ifu):
    global fncube, fnorigcube, fnmap, fncatalog

    dataroot = "../data"
    fncube="{dataroot}/sf2outcube_{field}_{ifu}.fits.gz".format(dataroot=dataroot,field=field,ifu=ifu)
    fnorigcube="{dataroot}/outcube_{field}_{ifu}.fits.gz".format(dataroot=dataroot,field=field,ifu=ifu)
    fnmap="{dataroot}/mmap_{field}_{ifu}.fits.gz".format(dataroot=dataroot,field=field,ifu=ifu)
    fncatalog="{dataroot}/msf2outcube_{field}_allifu.cat".format(dataroot=dataroot,field=field)

    return fncube, fnorigcube, fnmap, fncatalog

    
def load_data(field, ifu, vmax, fncatalog, fncube, fnorigcube, fnmap):
    global gbimages, gbimages_title, ii, current_record_idx , catalog, records, phz_pdz, phz_zz, c, oc, m


    h,t = os.path.split(fncatalog)
    newfncatalog = os.path.join( h , "m" + t)

    if os.path.exists(newfncatalog):
        print("Reading {}".format(newfncatalog))
        catalog = ascii.read(newfncatalog)
    else:
        print("Reading {}".format(fncatalog))
        catalog = ascii.read(fncatalog)


    B = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.B.original_psf.v2.fits".format(field=field,ifu=ifu)
    V = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.V.original_psf.v2.fits".format(field=field,ifu=ifu)
    gp = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.gp.original_psf.v2.fits".format(field=field,ifu=ifu)
    rp = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.rp.original_psf.v2.fits".format(field=field,ifu=ifu)
    ip = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.ip.original_psf.v2.fits".format(field=field,ifu=ifu)
    zp = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.zp.original_psf.v2.fits".format(field=field,ifu=ifu)
    Ks = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.Ks.original_psf.v5.fits".format(field=field,ifu=ifu)
    K = "../aux/thumbnails/{field}/ifu{ifu}_COSMOS.K.UV_original_psf.v1.fits".format(field=field,ifu=ifu)

    gbimages=[B,V,gp,rp,ip,zp,Ks,K]
    gbimages_title=["SUBARU B","SUBARU V","SUBARU gp","SUBARU rp","SUBARU ip","SUBARU zp","UltraVISTA Ks","CFHT K"]

    ii  = catalog["ifu"] == ifu
    ii *= ~np.isnan( catalog["manualx"] )
    ii *= catalog["N"] > 3
    current_record_idx = 0
    if not "class2" in catalog.colnames:
        #catalog.add_column(Column(["NA"]*len(catalog), name='class2', dtype='S100'), masked=False)
        catalog.add_column(Column(["NA"]*len(catalog), name='class2', dtype='S100'))
    if not "comment2" in catalog.colnames:
        #catalog.add_column(Column([""]*len(catalog), name='comment2', dtype='S500', masked=False))
        catalog.add_column(Column([""]*len(catalog), name='comment2', dtype='S500'))
    else:
        old_entries = list( catalog["comment2"] )
        catalog.remove_column("comment2")
        catalog.add_column( Column(old_entries, name='comment2', dtype='S500') )

    catalog = catalog.filled()
    records = catalog[ii]
    
    phz_pdz = fits.open("../pdz_cosmos2015_v1.3.fits.gz")
    c = spectrum.readSpectrum(fncube)
    oc = spectrum.readSpectrum(fnorigcube)
    m = spectrum.readSpectrum(fnmap)

    
    phz_zz = []
    for j in range(600):
         phz_zz.append( float( phz_pdz[1].header["TTYPE{}".format(j+75)][1:].replace("_",".") ) )
    ##record = records[0]
    #str(record["comment2"])  
    return records, gbimages
    
fits_cache = {}
def load_fits(fn):
    global cube_cache
    if not fn in fits_cache:
        hdulist = fits.open(fn)
        fits_cache[fn] = hdulist[0].header, hdulist[0].data
        hdulist.close()
    return fits_cache[fn]



def register_ds9staircase():
    # register color map
    from matplotlib.cm import register_cmap, cmap_d

    colors = []
    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*.3,kk*.3,kk*1)  )

    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*.3,kk*1,kk*.3)  )
    for ii in range(1,6):
        kk = ii/5.
        colors.append( (kk*1,kk*.3,kk*.3)  )
    colors = np.array(colors)
    xx = np.arange(len(colors), dtype=float)
    xx = xx/xx.max()

    ds9staircase = {'red': lambda v : np.interp(v, xx, colors[:,0]),
               'green': lambda v : np.interp(v, xx, colors[:,1]),
               'blue': lambda v : np.interp(v, xx, colors[:,2])}


    # Register all other colormaps
    register_cmap('ds9staircase', data=ds9staircase)

register_ds9staircase()



if not USEDS9:
    def updateCubes(fncube, fnorigcube, fnmap, vmax, cx, cy, cz, width=100):
        global figs
        cmap = cm.gray

        def showCube(fn, figN, vmin, vmax, cmap):
            plt.ioff()

            #hdulist = fits.open(fn)
            # we need this such that the images can transform to the same x,y coordinates
            
            header, data = load_fits(fn)
            cube_wcs = wcs.WCS(header)
            im = data[int(np.round(cz))]

            fig = plt.figure()

            ax = plt.subplot()            
            ax.imshow(im, origin = 'bottom', cmap = cmap, vmin=vmin, vmax=vmax)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.axvline(cx,color='w')
            ax.axhline(cy,color='w')
            ax.set_xlim([cx-width/2*imgscale,cx+width/2*imgscale])
            ax.set_ylim([cy-width/2*imgscale,cy+width/2*imgscale])
            
            with figN:
                clear_output(wait=True)
                display(ax.figure)
                
            return cube_wcs
                
        cube_wcs = showCube(fncube, figs[0], vmin=0., vmax=vmax, cmap=plt.get_cmap('ds9staircase')) 
        __ = showCube(fnorigcube, figs[1], vmin=0., vmax=vmax*3., cmap=plt.get_cmap('ds9staircase'))
        __ = showCube(fnmap, figs[2], vmin=0., vmax=1, cmap=cm.gray)
        
        return cube_wcs

                
    def plot_phz_ident(ax, rmatch = 2.5):
        global figs, phz_pdz, ii, catalog, current_record_idx
        
        records = catalog[ii]
        record = records[current_record_idx]
    
        ra,dec = record["ra_com"], record["dec_com"]

        dd = np.sqrt(((phz_pdz[1].data["RA"] - ra)*np.cos(np.deg2rad(dec)))**2. + (phz_pdz[1].data["DEC"] - dec)**2.)*3600.
        ii = dd < rmatch
        N = np.sum(ii)
        
        ax.plot([0], [0], 'bo', markersize=10)

        #print("{} objects within {} arcsec.".format( N, rmatch ) )

            
            
    def updateFigs(gbimages, cx, cy, cz, cube_wcs, width=100):
        global figs
        cmap = cm.gray_r
        
        for i,(fn,title) in enumerate(zip(gbimages, gbimages_title)):
            #print("Loading " + fn)

            header, data = load_fits(fn)
            im_wcs = wcs.WCS(header)
            im = data
            
            # translate cube x,y to ra/dec
            a,d,__ = [float(val) for val in cube_wcs.all_pix2world(cx,cy,cz,0)]

            imx, imy = [float(val) for val in im_wcs.all_world2pix(a,d,0)]

            plt.ioff()
            fig = plt.figure()

            # Create interval object
            #interval = MinMaxInterval()
            interval = PercentileInterval(97.)
            vmin, vmax = interval.get_limits(im)

            # Create an ImageNormalize object using a SqrtStretch object
            norm = ImageNormalize(vmin=vmin, vmax=vmax)
            
            plt.title(title)
            ax = plt.subplot()            
            ax.imshow(im, origin = 'bottom', norm=norm, cmap = cmap)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_xlim([imx-width/2*imgscale,imx+width/2*imgscale])
            ax.set_ylim([imy-width/2*imgscale,imy+width/2*imgscale])
            
            
            ax.axvline(imx,color='w')
            ax.axhline(imy,color='w')
            
            if i == 0:
                plot_phz_ident(ax)
            
            with figs[i+3]:
                clear_output(wait=True)
                display(ax.figure)

        
        
    def startup(field, ifu, fncube, fnorigcube, fnmap, gbimages, vmax):
        pass
        

    def clear_contours_on_images():
        print("Currently not available in non-ds9 version.")
        pass
    
    def transfer_contours():
        print("Currently not available in non-ds9 version.")
        pass
    
    def show_object(x,y,z,vmax,gbimages):
        plt.close('all')
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)
        updateFigs(gbimages, x, y, z, cube_wcs)
     

    def on_scale_up(b):
        global vmax
        vmax = vmax*1.05
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)

    def on_scale_down(b):
        global vmax
        vmax = vmax*.95
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)


    def on_prevslice_clicked(b):
        global x, y, z
        z = z  - 1.
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)

    def on_nextslice_clicked(b):
        global x, y, z
        z = z + 1.
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)
        
        
    def on_zoomin_clicked(b):
        global imgscale
        imgscale = imgscale*.66
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)
        updateFigs(gbimages, x, y, z, cube_wcs)
        
    def on_zoomout_clicked(b):
        global imgscale
        imgscale = imgscale*1.50
        cube_wcs = updateCubes(fncube, fnorigcube, fnmap, vmax, x, y, z)
        updateFigs(gbimages, x, y, z, cube_wcs)

        
  