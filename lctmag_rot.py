"""
LCTMag differential rotation — central meridian vφ profile
===========================================================
  Cell 1  imports + style
  Cell 2  paths / config
  Cell 3  I/O helpers
  Cell 4  fit helpers
  Cell 5  DATA LOADING   ← run once
  Cell 6  PLOT
"""


# %%  ── CELL 1: imports + style ────────────────────────────────────────────────

import os, glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import datetime, timedelta
from astropy.io import fits
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({
    "figure.facecolor": "white",  "axes.facecolor":  "white",
    "axes.edgecolor":   "#333",   "axes.labelcolor": "#222",
    "axes.linewidth":   0.8,      "axes.grid":       True,
    "grid.color":       "#e8e8e8","grid.linewidth":  0.5,
    "grid.linestyle":   "--",     "xtick.direction": "out",
    "ytick.direction":  "out",    "xtick.major.size": 4,
    "ytick.major.size": 4,        "xtick.color":     "#444",
    "ytick.color":      "#444",   "xtick.labelsize": 11,
    "ytick.labelsize":  11,       "font.family":     "sans-serif",
    "font.size":        11,       "axes.labelsize":  12,
    "axes.titlesize":   13,       "legend.fontsize": 10,
    "legend.framealpha": 0.9,     "legend.edgecolor": "#ccc",
    "figure.dpi":       150,      "savefig.dpi":     300,
    "savefig.bbox":     "tight",  "pdf.fonttype":    42,
})


# %%  ── CELL 2: paths / config ───────────────────────────────────────────────

BASE = ('/data/seismo/joshin/pipeline-test/local_correlation_tracking/'
        'pmi/diff_rot/data/')

DIRS = {
    "mag_4k": BASE + "data_for_im/data_mag_4k_2017/",
    "mag_2k": BASE + "data_for_im/data_mag_2k_2017/",
}

FILE_PATTERN = "*.hdf5"
VPHI_KEY     = "uphi"
LAT_KEY      = "latitude"
LON_KEY      = "longitude"

LAT_CLIP     = 90.0
LON_CLIP     = 90.0
LAT_FIT_MIN  = -90.0
LAT_FIT_MAX  =  90.0

SAVE_PREFIX  = "diffrot_mag"   # set None to display only


# %%  ── CELL 3: I/O helpers ────────────────────────────────────────────────────

def load_and_reduce(dataset_key, n_days=None, verbose=True):
    """
    Load HDF5 files → (n_days, n_lat, n_lon) cube → time-median map
    → central-meridian profile.

    Returns: lat, lon, vphi_map (n_lat × n_lon), vphi_prof (n_lat,)
    """
    files = sorted(glob.glob(os.path.join(DIRS[dataset_key], FILE_PATTERN)))
    if not files:
        raise FileNotFoundError(f"No files found for '{dataset_key}'")
    if n_days is not None:
        files = files[:n_days]

    with h5py.File(files[0], "r") as f:
        lat_full = f[LAT_KEY][:].astype(float)
        lon_full = f[LON_KEY][:].astype(float)

    lat_mask = np.abs(lat_full) <= LAT_CLIP
    lon_mask = np.abs(lon_full) <= LON_CLIP
    lat      = lat_full[lat_mask]
    lon      = lon_full[lon_mask]

    cube = np.full((len(files), lat.size, lon.size), np.nan)

    for i, fpath in enumerate(files):
        if verbose and i % 50 == 0:
            print(f"  [{dataset_key}] day {i+1}/{len(files)}")
        try:
            with h5py.File(fpath, "r") as f:
                vp = f[VPHI_KEY][:].astype(float)[:, lat_mask, :][:, :, lon_mask]
            cube[i] = np.nanmedian(vp, axis=0)   # collapse 4 intra-day frames
        except Exception as e:
            if verbose:
                print(f"  WARNING: skipping {fpath}: {e}")

    vphi_map  = np.nanmedian(cube, axis=0)                 # time median
    vphi_prof = vphi_map[:, int(np.argmin(np.abs(lon)))]  # central meridian

    return lat, lon, vphi_map, vphi_prof


# %%  ── CELL 4: fit helpers ─────────────────────────────────────────────────

def dr_model(lat_deg, A, B, C):
    """Ω(θ) = A + B sin²θ + C sin⁴θ"""
    s = np.sin(np.deg2rad(lat_deg))
    return A + B*s**2 + C*s**4


def fit_dr(lat, prof):
    mask = (lat >= LAT_FIT_MIN) & (lat <= LAT_FIT_MAX) & np.isfinite(prof)
    popt, pcov = curve_fit(dr_model, lat[mask], prof[mask],
                           p0=[2000., -300., -100.], maxfev=10_000)
    return popt, np.sqrt(np.diag(pcov))


# %%  ── CELL 5: DATA LOADING ───────────────────────────────────────────────────

print("Loading mag_4k ...")
lat_4k, lon_4k, map_4k, prof_4k = load_and_reduce("mag_4k", verbose=True)

print("\nLoading mag_2k ...")
lat_2k, lon_2k, map_2k, prof_2k = load_and_reduce("mag_2k", verbose=True)

popt_4k, perr_4k = fit_dr(lat_4k, prof_4k)
popt_2k, perr_2k = fit_dr(lat_2k, prof_2k)

print(f"\nmag_4k  A={popt_4k[0]:.1f}  B={popt_4k[1]:.1f}  C={popt_4k[2]:.1f}")
print(f"mag_2k  A={popt_2k[0]:.1f}  B={popt_2k[1]:.1f}  C={popt_2k[2]:.1f}")


# %%  ── CELL 6: PLOT ────────────────────────────────────────────────────────────
#
#  ── things you might want to change quickly ───────────────────────────────────
#
#     DATA / PROFILES     lat_4k, prof_4k, popt_4k  (loaded above)
#     COLORS              c_4k, c_2k
#     AXIS LIMITS         ax.set_xlim / ax.set_ylim
#     EXTRA LINES         e.g. ax.axhline(...) for a reference rotation rate
#     FIT CURVE           lat_fit, dr_model(lat_fit, *popt_4k)
#     LEGEND / LABELS     ax.set_xlabel, ax.legend, ax.set_title
#
# ─────────────────────────────────────────────────────────────────────────────

c_4k = "#c0392b"    # HMI  4K  colour
c_2k = "#e8a090"    # PMI  2K  colour



rot_ref = 456.03128739456554 # [nHz]

class Rot(object):
    def __init__(self, infile):
        hdus = fits.open(infile)
        self.nHz_hemi = hdus['rot_nHz'].data
        nt, nr, nlat = self.nHz_hemi.shape
        # coordinates
        self.lat_hemi = hdus['latitude'].data['lat']
        self.r = hdus['rmesh'].data['r']
        self.tstart = np.array([datetime.strptime(s, '%Y-%m-%d')
                for s in hdus['time'].data['tstart']])
        self.dats = self.tstart + timedelta(days=72/2)
        # pad
        self.nHz = np.pad(self.nHz_hemi, ((0,0),(0,0),(0,nlat-1)), mode='reflect') - rot_ref
        self.lat = np.pad(self.lat_hemi, (0,nlat-1), mode='reflect')
        self.lat[:nlat] *= -1
        self.ms = self.nHz*1e-9*2*np.pi*696e6*np.cos(np.deg2rad(self.lat[None,None,:]))


indir = Path('/data/seismo/zhichao/codes/proj_drot/dwnld-globalhs/output')
period = '2010.04.30-2025.09.14'
infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF36_{period}.fits'
rot_36 = Rot(infile)
sinlat = np.sin(np.deg2rad(rot_36.lat))
rsinth = 696e6*np.cos(np.deg2rad(rot_36.lat))
rot_ms_36 = np.nanmean(rot_36.ms, axis=0)


infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF18_{period}.fits'
rot_18 = Rot(infile)
sinlat = np.sin(np.deg2rad(rot_18.lat))
rsinth = 696e6*np.cos(np.deg2rad(rot_18.lat))
rot_ms_18 = np.nanmean(rot_18.ms, axis=0)

infile = indir / f'hmi.v_sht_2drls-72day-NACOEFF6_{period}.fits'
rot_6 = Rot(infile)
sinlat = np.sin(np.deg2rad(rot_6.lat))
rsinth = 696e6*np.cos(np.deg2rad(rot_6.lat))
rot_ms_6 = np.nanmean(rot_6.ms, axis=0)

r_targ = [0.97]


lat_fit = np.linspace(LAT_FIT_MIN, LAT_FIT_MAX, 500)

fig, ax = plt.subplots(figsize=(6, 5))

# ── data (translucent) ────────────────────────────────────────────────────────
ax.plot(lat_4k, prof_4k, color=c_4k, lw=0, marker="o",
        ms=2.5, zorder=1)
ax.plot(lat_2k, prof_2k, color=c_2k, lw=0, marker="o",
        ms=2.5, zorder=1)

# --- global helioseismology (solid) ────────────────────────────────────────────────────────
for k, r1 in enumerate(r_targ):
    ir_36 = np.argmin(np.abs(rot_36.r - r1))
    ir_18 = np.argmin(np.abs(rot_18.r - r1))
    ir_6  = np.argmin(np.abs(rot_6.r - r1))
    ax.plot(rot_36.lat, rot_ms_36[ir_36], color="#2980b9", lw=1, label=f"Global helioseismology (NA={36})")
    ax.plot(rot_18.lat, rot_ms_18[ir_18], color="#3498db", lw=1, label=f"Global helioseismology (NA={18})")
    ax.plot(rot_6.lat,  rot_ms_6[ir_6],  color="#5dade2", lw=1, label=f"Global helioseismology (NA={6})")



# ── reference line ────────────────────────────────────────────────────────────
ax.axhline(0, color="#aaaaaa", lw=0.8, ls=":")

# ── labels / limits ───────────────────────────────────────────────────────────
ax.set_xlabel("Latitude (°)")
ax.set_ylabel(r"$v_\varphi\ \mathrm{(m\,s^{-1})}$")
ax.set_title(r"LCTMag — differential rotation, central meridian")
ax.set_xlim(LAT_FIT_MIN, LAT_FIT_MAX)
ax.xaxis.set_major_locator(plt.MultipleLocator(15))
ax.legend(fontsize=8)

fig.tight_layout()

# if SAVE_PREFIX:
#     fig.savefig(f"{SAVE_PREFIX}.pdf")
#     print(f"Saved → {SAVE_PREFIX}.pdf")

plt.show()

# %%  ── CELL 6: PLOT ────────────────────────────────────────────────────────────
#
#  Units: nHz  (Ω - Ω_ref), as is standard in helioseismology.
#
#  Conversion for LCT profiles (given in m/s):
#      Ω [nHz] = v_φ / (2π · R_☉ · cos λ) × 1e9
#  where λ is latitude in radians and R_☉ = 696 Mm.
#  The helio data is already in nHz (rot_obj.nHz), so we use that directly
#  rather than the pre-converted .ms arrays.
#
#  ── quick-change handles ──────────────────────────────────────────────────────
R_TARG  = [1.00, 0.98, 0.96]
NCOLS   = 3

C_4K = "#117a8b"   # deep teal — HMI 4K
C_2K = "#76d7ea"   # light teal — PMI 2K
C_NA = {36: "#b7440a", 18: "#e67e22", 6: "#f5cba7"}
LW_NA   = 1.4

XLIM    = (LAT_FIT_MIN, LAT_FIT_MAX)
R_SUN   = 696.0e6          # metres
ROT_REF = 0  # nHz  (same reference subtracted in Rot class)
# ─────────────────────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D


def vphis_to_nHz(lat_deg, vphi_ms):
    """
    Convert v_φ [m/s] → differential rotation rate [nHz] relative to ROT_REF.

        Ω_total [nHz] = v_φ / (2π · R_☉ · cos λ) × 1e9
        Ω_diff  [nHz] = Ω_total - ROT_REF

    Poles (|λ| → 90°) are masked to avoid division by near-zero cos λ.
    """
    cos_lat = np.cos(np.deg2rad(lat_deg))
    with np.errstate(invalid="ignore", divide="ignore"):
        omega_total = vphi_ms / (2.0 * np.pi * R_SUN * cos_lat) * 1e9
    omega_diff = omega_total - ROT_REF
    # mask latitudes where cos λ < 0.1 (|λ| > ~84°) — numerically unreliable
    omega_diff[np.abs(lat_deg) > 85] = np.nan
    return omega_diff


# ── convert LCT profiles once ─────────────────────────────────────────────────
prof_4k_nHz = vphis_to_nHz(lat_4k, prof_4k)
prof_2k_nHz = vphis_to_nHz(lat_2k, prof_2k)

# ── helio: time-mean nHz at each radius (Ω - ROT_REF already in rot_obj.nHz) ─
rot_nHz_mean = {
    36: np.nanmean(rot_36.nHz, axis=0),   # shape (nr, nlat_full)
    18: np.nanmean(rot_18.nHz, axis=0),
    6:  np.nanmean(rot_6.nHz,  axis=0),
}
rot_objs = {36: rot_36, 18: rot_18, 6: rot_6}

# ─────────────────────────────────────────────────────────────────────────────
lat_fit  = np.linspace(*XLIM, 500)
nrows    = int(np.ceil(len(R_TARG) / NCOLS))
npanels  = len(R_TARG)

fig = plt.figure(figsize=(NCOLS * 4.2, nrows * 3.8))
gs  = gridspec.GridSpec(nrows, NCOLS, figure=fig,
                        hspace=0.42, wspace=0.32)

for k, r1 in enumerate(R_TARG):
    row, col = divmod(k, NCOLS)
    ax = fig.add_subplot(gs[row, col])

    # ── global helioseismology ────────────────────────────────────────────────
    for na in [36, 18, 6]:
        rot_obj = rot_objs[na]
        ir      = np.argmin(np.abs(rot_obj.r - r1))
        # nHz already relative to ROT_REF; select surface hemisphere then mirror
        ax.plot(rot_obj.lat, rot_nHz_mean[na][ir],
                color=C_NA[na], lw=LW_NA, zorder=2,
                label=rf"$N_A = {na}$")

    # ── LCTMag profiles ───────────────────────────────────────────────────────
    ax.plot(lat_4k, prof_4k_nHz,
            color=C_4K, lw=0, marker="o", ms=2.2, alpha=0.55, zorder=3)
    ax.plot(lat_2k, prof_2k_nHz,
            color=C_2K, lw=0, marker="o", ms=2.2, alpha=0.55, zorder=3)

    # ── panel label ───────────────────────────────────────────────────────────
    ax.text(0.04, 0.96, rf"$r = {r1:.2f}\,R_\odot$",
            transform=ax.transAxes, fontsize=12,
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", fc="white",
                      ec="#cccccc", alpha=0.85))

    ax.axhline(0, color="#bbbbbb", lw=0.7, ls=":")
    ax.set_xlim(*XLIM)
    ax.set_ylim(-200, 20)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.xaxis.set_major_locator(mticker.MultipleLocator(30))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.tick_params(labelsize=12)

    if col == 0:
        ax.set_ylabel(r"$\Omega - \Omega_{\rm ref}\ \mathrm{(nHz)}$", fontsize=12)
    if row == nrows - 1 or k == npanels - 1:
        ax.set_xlabel("Latitude (°)", fontsize=12)

    ax.legend(fontsize=10, loc="lower center",
              handlelength=1.4, handletextpad=0.5,
              borderpad=0.5, labelspacing=0.35)

# ── hide unused panels ────────────────────────────────────────────────────────
for k in range(npanels, nrows * NCOLS):
    fig.add_subplot(gs[divmod(k, NCOLS)]).set_visible(False)

# ── shared LCT legend (bottom) ────────────────────────────────────────────────
lct_handles = [
    Line2D([0], [0], color=C_4K, lw=0, marker="o", ms=5,
           label="LCTMag HMI 4K"),
    Line2D([0], [0], color=C_2K, lw=0, marker="o", ms=5,
           label="LCTMag PMI 2K"),
]
fig.legend(handles=lct_handles,
           loc="lower center", ncol=2, fontsize=9,
           framealpha=0.9, edgecolor="#cccccc",
           bbox_to_anchor=(0.5, -0.14))

fig.suptitle(
    r"Differential rotation: LCTMag vs global helioseismology"
    r"$\quad(\Omega_{\rm ref} = 456.03\ \mathrm{nHz})$",
    fontsize=16, y=1.01,
)

fig.savefig(f"{SAVE_PREFIX}_multiradius_nHz.pdf")
# print(f"Saved → {SAVE_PREFIX}_multiradius_nHz.pdf")

plt.show()
# %%
