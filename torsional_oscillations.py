# %%
import os
from datetime import datetime, timedelta
from pathlib import Path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
from matplotlib.dates import date2num, DateFormatter, YearLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec

class Rot():
     NU_CR = 1e9/(25.38*86400) # 456.03128739456554 [nHz]

     def __init__(self, infile):
         hdus = fits.open(infile)
         self.nHz_hemi = hdus['rot_nHz'].data
         nt, nr, nlat = self.nHz_hemi.shape
         print('read', infile)
         self.fname = os.path.basename(infile)
         # coordinates
         self.lat_hemi = hdus['latitude'].data['lat']
         self.r = hdus['rmesh'].data['r']
         self.tstart = np.array([datetime.strptime(s, '%Y-%m-%d')
                 for s in hdus['time'].data['tstart']])
         self.dats = self.tstart + timedelta(days=72/2)
         self.xdat = date2num(self.dats)
         # pad
         self.nHz = np.pad(self.nHz_hemi, ((0,0),(0,0),(0,nlat-1)), mode='reflect') - Rot.NU_CR
         self.lat = np.pad(self.lat_hemi, (0,nlat-1), mode='reflect')
         self.lat[:nlat] *= -1
         self.ms = self.nHz*1e-9*2*np.pi*696e6*np.cos(np.deg2rad(self.lat[None,None,:]))
         # extent
         xl, xh = self.xdat[0]-72/2, self.xdat[-1]+72/2
         dlat = self.lat[1] - self.lat[0]
         yl, yh = self.lat[0]-0.5*dlat, self.lat[-1]+0.5*dlat
         self.extent = [xl, xh, yl, yh]

def rmtmean(data, xdat, xlim): # data[lat,t]
     trng = (xdat >= xlim[0]) & (xdat <= xlim[1])
     return data - np.nanmean(data[:,trng], axis=1)[:,None]

indir = Path('/data/seismo/zhichao/codes/proj_drot/dwnld-globalhs/output')

period = '2010.04.30-2025.09.14'
infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF36_{period}.fits'
na36 = Rot(infile)
infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF18_{period}.fits'
na18 = Rot(infile)
infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF6_{period}.fits'
na6 = Rot(infile)

# %%
## plot
nrows, ncols = 3, 1
gs = GridSpec(nrows, ncols, hspace=0.5)
fig = plt.figure(figsize=(6,8))
fig.suptitle(f'$u_\\phi - \\langle u_\\phi [file://\\rangle$]\\rangle$', y=0.95, fontsize='xx-large')
for irow, (rot,r1,clim) in enumerate([
     (na36, 0.99, [-10,10]),
     (na18, 0.99, [-10,10]),
     (na6,  0.99, [-10,10]),
     # (na36, 0.95, [-10,10]),
     # (na18, 0.95, [-10,10]),
     # (na6,  0.95, [-10,10]),
     ]):
     # rot.ms[t,r,lat]
     # ux1[lat,t]
     ir = np.argmin(abs(rot.r - r1))
     ux1 = np.swapaxes(rot.ms[:,ir,:], 0, 1)
     ux1 = rmtmean(ux1, rot.xdat, rot.extent[:2])

     ax = fig.add_subplot(gs[irow])
     im = ax.imshow(ux1, origin='lower', interpolation='nearest', cmap='bwr', clim=clim, aspect='auto', extent=rot.extent)
     ax.axhline(c='k', lw=0.7)
     # colorbar
     l, b, w, h = ax.get_position().bounds
     cax = fig.add_axes([l+1.02*w, b, 0.02*w, h])
     fig.colorbar(im, cax=cax, ticks=clim)
     # cax.set_ylabel('$u_\\phi - \\langle u_\\phi [file://\\rangle$]\\rangle$', labelpad=-18, fontsize='x-large')
     cax.set_title('[m/s]', y=1.05)
     # decoration
     ax.xaxis.set_major_locator(YearLocator(5))
     ax.xaxis.set_minor_locator(YearLocator(1))
     ax.yaxis.set_major_locator(MultipleLocator(30))
     ax.yaxis.set_minor_locator(MultipleLocator(10))
     ax.xaxis.set_major_formatter(DateFormatter('%Y'))
     yticks = ax.get_yticks()
     ax.set_yticks(yticks)
     ax.set_yticklabels([f'${v:g}\\degree$' for v in yticks])
     ax.set_xlim(rot.extent[:2])
     ax.set_ylim(-90,90)
     ax.set_xlabel('Year')
     ax.set_ylabel('Latitude')
     i1 = rot.fname.find('NACOEFF')
     i2 = rot.fname[i1:].find('_')
     ax.set_title(f'{rot.fname[i1:i1+i2]}, $r/R_\\odot = {r1:.2f}$')