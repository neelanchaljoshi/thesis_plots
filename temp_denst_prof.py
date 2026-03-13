# %%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation, AutoMinorLocator
from scipy.constants import k, m_u

# --- Load VAL3C data ---
height, T, pe_arr, pg = np.loadtxt("VAL3C_rutten.dat",
    unpack=True,
    skiprows=1,
    usecols=(0, 3, 4, 5))

pg     /= 10
pe_arr /= 10

# --- Density ---
mu      = 1.0 / (1.0 + pe_arr / pg)
rho     = (mu * m_u * pg) / (k * T)
rho_cgs = rho * 1e-3

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "serif",
    "font.serif":         ["Palatino Linotype", "Palatino", "Book Antiqua",
                           "Georgia", "DejaVu Serif"],
    "mathtext.fontset":   "cm",
    "axes.linewidth":     1.1,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "xtick.major.size":   6,
    "ytick.major.size":   6,
    "xtick.minor.size":   3.5,
    "ytick.minor.size":   3.5,
    "xtick.top":          False,
    "ytick.right":        False,
    "figure.dpi":         150,
})

# Region colours — all distinct, no blues
C_photo  = "#2a7d4f"   # deep malachite green
C_chrom  = "#b5480f"   # dark burnt orange
C_trans  = "#7b2d8b"   # rich purple
C_corona = "#a6122d"   # deep crimson

# Data line colours
C_temp   = "#1464a0"   # cerulean blue
C_rho    = "#c47a00"   # dark golden amber

Tmin = 3_000
Tmax = 35_000

fig, ax = plt.subplots(figsize=(9, 5.5))
fig.subplots_adjust(left=0.09, right=0.87, top=0.93, bottom=0.13)

# Restore top spine explicitly
ax.spines["top"].set_visible(True)

# ── Region shading + labels (no dark edge bands) ──────────────────────────────
regions = [
    (height.min(), 500,          C_photo,  "Photosphere",        False),
    (500,          2100,         C_chrom,  "Chromosphere",       False),
    (2100,         2300,         C_trans,  "Transition\nRegion", True),
    (2300,         height.max(), C_corona, "Corona",             True),
]
for x0, x1, color, label, vert in regions:
    ax.axvspan(x0, x1, alpha=0.18, color=color, lw=0, zorder=0)
    if x0 != height.min():
        ax.axvline(x=x0, color=color, lw=1.0, ls=(0, (5, 3)), alpha=0.7, zorder=2)
    xmid = (x0 + x1) / 2
    ax.text(xmid, 24_500, label,
            ha="center", va="center",
            fontsize=8, color=color,
            rotation=90 if vert else 0,
            fontweight="bold", alpha=0.95)

# ── Temperature ───────────────────────────────────────────────────────────────
ax.plot(height, T, color=C_temp, lw=2.2, solid_capstyle="round", zorder=5)

ax.set_xlabel("Height [km]", fontsize=11, labelpad=6)
ax.set_ylabel("Temperature [K]", fontsize=11, color=C_temp, labelpad=6)
ax.set_yscale("log")
ax.set_xlim([height.min(), height.max()])
ax.set_ylim([Tmin, Tmax])
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(
    lambda x, _: rf"$\mathregular{{{int(x):,}}}$"))
ax.yaxis.set_minor_formatter(mpl.ticker.FuncFormatter(
    lambda x, _: rf"$\mathregular{{{int(x):,}}}$"
    if x in (5000, 10000, 20000) else ""))

# Set tick colours after formatters so labels aren't recreated with default colour
ax.tick_params(axis="y", labelcolor=C_temp, labelsize=9, which="both")
ax.tick_params(axis="x", labelsize=9)
plt.setp(ax.get_yticklabels(which="both"), color=C_temp)

# ── Density ───────────────────────────────────────────────────────────────────
ax2 = ax.twinx()
ax2.spines["top"].set_visible(True)
ax2.plot(height, rho_cgs,
         color=C_rho, lw=2.2, ls=(0, (6, 2)),
         solid_capstyle="round", zorder=5)
ax2.set_ylabel(r"Density  [g cm$^{-3}$]", fontsize=11, color=C_rho, labelpad=6)
ax2.tick_params(axis="y", labelcolor=C_rho, direction="out", labelsize=9)
ax2.set_yscale("log")
ax2.set_ylim([rho_cgs.min() * 0.4, rho_cgs.max() * 4])
ax2.yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))

# ── Legend ────────────────────────────────────────────────────────────────────
leg_handles = [
    mpl.lines.Line2D([0],[0], color=C_temp, lw=2.2,           label="Temperature"),
    mpl.lines.Line2D([0],[0], color=C_rho,  lw=2.2, ls=(0,(6,2)), label="Density"),
]
ax.legend(handles=leg_handles,
          loc="lower center", bbox_to_anchor=(0.5, 0.02),
          fontsize=9.5, framealpha=0.95,
          edgecolor="#aaaaaa", handlelength=2.8,
          borderpad=0.7)

plt.savefig("temp_dens_profile.pdf", bbox_inches="tight")
# plt.savefig("temp_dens_profile.png", dpi=300, bbox_inches="tight")
plt.show()
# %%

# %%
