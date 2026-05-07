#!/usr/bin/env python3
# %%
from pathlib import Path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 24
from matplotlib.patches import Arc

def imshow_contour_stack(ax, X0, Y0, dx, Nt, Xspace, dx2, imgs, contours,
        lev=None, cticks=None, cntrlev=None,
        titles=[], clab='', unit='', xoff=0, yoff=0, cmap='seismic',
        ix_label_r=0, ix_label_lat=-1, ix_colorbar=-1):
    cmin, cmax = (np.nan, np.nan) if lev is None else lev[[0,-1]]
    zmin, zmax = np.nanmin(imgs), np.nanmax(imgs)
    if zmin < cmin and zmax > cmax:
        extend = 'both'
    elif zmin < cmin and zmax < cmax:
        extend = 'min'
    elif zmin > cmin and zmax > cmax:
        extend = 'max'
    else:
        extend = 'neither'

    nx = len(imgs)
    for ix,(img,cntr) in enumerate(zip(imgs,contours)):
        if ix < Nt:
            DX = dx*ix
            DY = 0
            S = 1
        else:
            DX = dx*(Nt-1) + Xspace + (ix-Nt)*dx2/2.
            #DY = -0.5
            DY = 0
            S = 0.5
        X = X0*S + DX
        Y = Y0*S + DY
        im = ax.contourf(X, Y, img, cmap=cmap, levels=lev, extend=extend, zorder=1-10*ix)
        if cntrlev is not None:
            cs = ax.contour(X, Y, cntr, cntrlev, zorder=2-10*ix,
                    colors='k', linestyles='solid')

        # --- solar radii ---
        for arc_r in [1, 0.9, 0.8, 0.7]:
            ls = '-' if arc_r == 1 or arc_r == 0.7 else ':'
            ax.add_patch(Arc([DX,DY], 2*arc_r*S, 2*arc_r*S, angle=0, theta1=-90, theta2=90, color='k', ls=ls, zorder=2-10*ix))
            if ix == ix_label_r and (arc_r == 1 or arc_r == 0.7):
                ax.text(-0.05+DX, DY+arc_r*S, '%.1f' % (arc_r), fontsize='medium', ha='right', va='center')
                ax.text(-0.05+DX, DY-arc_r*S, '%.1f' % (arc_r), fontsize='medium', ha='right', va='center')
        if ix == ix_label_r:
            ax.text(-0.1, 0.5, '$r/\\mathrm{R}_\\odot$', ha='center', va='center', fontsize='large')
        if ix < Nt:
            ax.plot([DX,DX], [0.7,1], 'k', zorder=2-10*ix)
            ax.plot([DX,DX], [-0.7,-1], 'k', zorder=2-10*ix)
        else:
            ax.plot([DX,DX], [0.7*S+DY,1*S+DY], 'k', zorder=2-10*ix)
            ax.plot([DX,DX], [-0.7*S+DY,-1*S+DY], 'k', zorder=2-10*ix)

        # --- latitude ---
        for lat in range(-75,90,15):
            ang = np.deg2rad(lat)
            x2 = np.cos(ang)
            y2 = np.sin(ang)
            x1, y1 = 0.7*x2, 0.7*y2
            ax.plot([x1*S+DX,x2*S+DX], [y1*S+DY,y2*S+DY], ':k', zorder=2-10*ix)
            if ix == ix_label_lat:
                if abs(lat) >= 60: continue
                if lat > 0: hemi = 'N'
                elif lat < 0: hemi = 'S'
                else:hemi = ''
                ax.text(1.2*x2*S+DX, 1.2*y2*S+DY, '$%g\\degree$%s' % (abs(lat), hemi),
                        rotation=lat, va='center', ha='center', fontsize='medium',
                        bbox=dict(fc='none', ec='none'))

        # --- colorbar ---
        ax.set_aspect('equal')
        if ix == ix_colorbar:
            ax.apply_aspect()
            l, b, w, h = ax.get_position().bounds
            cax = fig.add_axes([l-1.25*w, b+0.1*h, 0.03*w, 0.8*h])
            cbar = fig.colorbar(im, cax=cax, orientation='vertical', ticks=cticks)
            cax.yaxis.set_ticks_position('left')
            cax.text(-7+xoff, 0.40+yoff, clab, fontsize='xx-large', rotation=90, va='top', ha='center')
            cax.text(-7+xoff, 0.44+yoff, unit, fontsize='large', rotation=90, va='bottom', ha='center')

        if len(titles) > 0:
            if ix < Nt:
                ax.text(DX, 1.1, titles[ix], fontsize='large')
            else:
                if ix % 2 == 0:
                    ax.text(DX+0.2, 0.55, titles[ix], ha='center', va='bottom')
                else:
                    ax.text(DX+0.2,-0.6, titles[ix], ha='center', va='top')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(left=False, right=False, top=False, bottom=False,
            which='both', labelleft=False, labelbottom=False)
    ax.set_xlim(-0.01, 1.01*S + DX)
    ax.set_ylim(-1.01, 1.01)

class Data(object):
    def __init__(self, infile, is_cyc=False):
        hdus = fits.open(infile)
        zeta_phi = hdus['zeta_phi'].data
        psi_scal = hdus['psi_scal'].data
        uth = hdus['uth'].data
        ur = hdus['ur'].data

        rsun = hdus['rmesh'].header['rsun'] # [m]
        rb = hdus['rmesh'].header['rb'] # [rsun]
        r = hdus['rmesh'].data['r'] # [rsun]
        th = hdus['theta'].data['th'] # [rad]

        if is_cyc:
            uth[r < rb,:] = np.nan
            ur[r < rb,:] = np.nan
            psi_scal[r < rb,:] = np.nan
            psi_scal[:,[0,-1]] = psi_scal[:,[1,-2]]
            nr, nth = uth.shape # (178, 181)
        else:
            uth[:,r < rb,:] = np.nan
            ur[:,r < rb,:] = np.nan
            psi_scal[:,r < rb,:] = np.nan
            psi_scal[:,:,[0,-1]] = psi_scal[:,:,[1,-2]]
            nt, nr, nth = uth.shape # (?, 178, 181)
            periods = hdus['time'].data['period']
            decimalyears = hdus['time'].data['decimalyear']
            yrs = np.asarray([int(v) for v in decimalyears])
            assert len(periods) == nt
            self.yrs = yrs
            self.decimalyears = decimalyears
            self.periods = periods
        assert len(r) == nr
        assert len(th) == nth
        self.uth = uth
        self.ur = ur
        self.psi = psi_scal
        self.r = r
        self.th = th

def new_zero(arr, new):
    l = list(arr)
    assert l.count(0) == 1
    arg = l.index(0)
    l.insert(arg+1, new)
    l.insert(arg, -new)
    l.remove(0)
    return np.asarray(l)

def scaled(psi):
    vmax = np.nanmax(np.abs(psi))
    return psi/vmax

npow = 2
def powered(arr):
    return np.sign(arr)*np.abs(arr)**npow
def rooted(arr):
    return np.sign(arr)*np.abs(arr)**(1./npow)

# %%
## getdata
# indir = '/scratch/seismo/zhichao/proj_mflow3/inverted-mflow,v5/daily_mdi+gong_ratio'
indir = Path('/data/seismo/zhichao/codes/proj_mflow3/python3/data')
infile = indir / 'psi_MDIGONG-daily_RegCurl_alpha65.fits'
mrg = Data(infile)

infile = indir / 'psi_MDI2003MayGONG-daily_cycle23_alpha70_RegCurl.fits'
cy23 = Data(infile, is_cyc=True)
infile = indir / 'psi_MDI2003MayGONG-daily_cycle24_alpha70_RegCurl.fits'
cy24 = Data(infile, is_cyc=True)

## remap
assert np.all(mrg.r == cy23.r)
assert np.all(mrg.th == cy23.th)
periods = mrg.periods
yrs = mrg.yrs
r = mrg.r
th = mrg.th
R = r[:,None]
X = R*np.sin(th)
Y = R*np.cos(th)

indx = list(range(0,len(yrs),3)) + [-1]
#titles = ['Cyc 23', 'Cyc 24'] + ['\n$-$'.join(s.split('-')) for s in periods[indx]]
titles = ['Cycle 23', 'Cycle 24'] + ['\n$\\downarrow$\n'.join(s.split('-')) for s in periods[indx]]
contour_param = (scaled(mrg.psi[indx]), scaled(cy23.psi), scaled(cy24.psi), powered(np.linspace(-1,1,9)))
params = [
        # quantity for each year, quantity for cycles 23 and 24,
        # lev, cticks, clab, unit
        (mrg.uth[indx], cy23.uth, cy24.uth,
            new_zero(np.linspace(-12,12,17),0.6), np.linspace(-12,12,9),
            '$u_\\theta$', ' (m s$^{-1}$)'),
        (scaled(mrg.psi[indx]), scaled(cy23.psi), scaled(cy24.psi),
            powered(np.linspace(-1,1,17)), powered(np.linspace(-1,1,5)), # powered(np.linspace(-1,1,9)),
            '$\\psi$', ' (normalized)'),
        ]

# %%
## plot
nrows = len(params)
fig = plt.figure(figsize=(20,16))

# --- top ---
q1, q23, q24, lev, cticks, clab, unit = params[0]
# panels A, C
ax = fig.add_subplot(2,1,1)
nt = 2
imgs = [q23, q24] + list(q1)
contours = [None, None] + list(q1)
imshow_contour_stack(ax, X, Y, 0.8, nt, 1.5, 0.8, imgs, contours,
        lev=lev, cticks=cticks, cntrlev=None,
        clab=clab, unit=unit, xoff=-3, yoff=-2,
        titles=titles, cmap='bwr_r',
        ix_label_r=0, ix_label_lat=1, ix_colorbar=1)
ax.text(-1.1, 1.1, 'A', fontsize='xx-large', weight='bold')
ax.text(2.2, 1.1, 'C', fontsize='xx-large', weight='bold')

# --- bottom ---
q1, q23, q24, lev, cticks, clab, unit = params[1]
q1_cntr, q23_cntr, q24_cntr, cntrlev = contour_param
# panels B, D
ax = fig.add_subplot(2,1,2)
nt = 2
imgs = [q23, q24] + list(q1)
contours = [q23_cntr, q24_cntr] + list(q1_cntr)
imshow_contour_stack(ax, X, Y, 0.8, nt, 1.5, 0.8, imgs, contours,
        lev=lev, cticks=cticks, cntrlev=cntrlev,
        clab=clab, unit=unit, xoff=-3, yoff=-0.7,
        titles=titles, cmap='coolwarm',
        ix_label_r=0, ix_label_lat=1, ix_colorbar=1)
ax.text(-1.1, 1.1, 'B', fontsize='xx-large', weight='bold')
ax.text(2.2, 1.1, 'D', fontsize='xx-large', weight='bold')

# %%
## plot
fig = plt.figure(figsize=(20,8))  # reduced height since only one row now

# --- top (only row now) ---
q1, q23, q24, lev, cticks, clab, unit = params[0]
ax = fig.add_subplot(1,1,1)
nt = 2
imgs = [q23, q24] + list(q1)
# contours = [None, None] + list(q1)
cntrlev = np.linspace(-12, 12, 9)   # or whatever levels you want
contours = [q23, q24] + list(q1)
imshow_contour_stack(ax, X, Y, 0.8, nt, 1.5, 0.8, imgs, contours,
        lev=lev, cticks=cticks, cntrlev=cntrlev,
        clab=clab, unit=unit, xoff=-3, yoff=-2,
        titles=titles, cmap='bwr',
        ix_label_r=0, ix_label_lat=1, ix_colorbar=1)
# fig.savefig('merid_circ_zcl.pdf', bbox_inches='tight')

# %%
