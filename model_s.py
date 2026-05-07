# # %%
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# from matplotlib.ticker import AutoMinorLocator, LogFormatterSciNotation

# # --- Load solar model ---
# url = "https://users-phys.au.dk/~jcd/solar_models/cptrho.l5bi.d.15c"
# data = np.genfromtxt(url, comments="#")

# radius = data[:,0]      # Radius (R/Rsun)
# rho_cgs = data[:,2]    # Density (g/cm^3)
# T = data[:,5]      # Temperature (K)

# # ── Style ────────────────────────────────────────────────────────────────────
# plt.rcParams.update({
#     "font.family":        "serif",
#     "font.serif":         ["Palatino Linotype", "Palatino", "Book Antiqua",
#                            "Georgia", "DejaVu Serif"],
#     "mathtext.fontset":   "cm",
#     "axes.linewidth":     1.1,
#     "xtick.direction":    "out",
#     "ytick.direction":    "out",
#     "xtick.major.size":   6,
#     "ytick.major.size":   6,
#     "xtick.minor.size":   3.5,
#     "ytick.minor.size":   3.5,
#     "xtick.top":          False,
#     "ytick.right":        False,
#     "figure.dpi":         150,
# })

# # Region colours — all distinct, no blues
# C_photo  = "#2a7d4f"   # deep malachite green
# C_chrom  = "#b5480f"   # dark burnt orange
# C_trans  = "#7b2d8b"   # rich purple
# C_corona = "#a6122d"   # deep crimson

# # Data line colours
# C_temp   = "#1464a0"   # cerulean blue
# C_rho    = "#c47a00"   # dark golden amber


# # --- Colors ---
# # C_temp  = "crimson"
# # C_rho   = "royalblue"
# C_core  = "#FFD700"  # gold
# C_rad   = "#FFA500"  # orange
# C_cz    = "#87CEEB"  # skyblue

# # --- Axis limits ---
# Tmin, Tmax = T.min()*0.8, T.max()*1.05
# rmin, rmax = radius.min(), radius.max()*1.01

# # --- Figure ---
# fig, ax = plt.subplots(figsize=(9, 5.5))
# fig.subplots_adjust(left=0.09, right=0.87, top=0.93, bottom=0.13)

# ax.spines["top"].set_visible(True)

# # ── Region shading ────────────────────────────────────────────────────────────
# regions = [
#     (radius.min(), 0.25, C_core,  "Core",            False),
#     (0.25,          0.71, C_rad,   "Radiative\nZone", False),
#     (0.71,          radius.max(), C_cz, "Convection\nZone", False),
# ]

# for x0, x1, color, label, vert in regions:
#     ax.axvspan(x0, x1, alpha=0.18, color=color, lw=0, zorder=0)

#     if x0 != radius.min():
#         ax.axvline(x=x0, color=color, lw=1.0, ls=(0,(5,3)), alpha=0.7, zorder=2)

#     xmid = (x0 + x1)/2
#     ax.text(xmid, Tmax*0.35, label,
#             ha="center", va="center",
#             fontsize=8, color=color,
#             rotation=90 if vert else 0,
#             fontweight="bold", alpha=0.95)

# # ── Temperature ───────────────────────────────────────────────────────────────
# ax.plot(radius, T, color=C_temp, lw=2.2, solid_capstyle="round", zorder=5)

# ax.set_xlabel(r"Radius [$R_\odot$]", fontsize=11, labelpad=6)
# ax.set_ylabel("Temperature [K]", fontsize=11, color=C_temp, labelpad=6)

# ax.set_yscale("log")
# ax.set_xlim([radius.min(), radius.max()])
# ax.set_ylim([Tmin, Tmax])

# ax.xaxis.set_minor_locator(AutoMinorLocator(5))

# ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(
#     lambda x, _: rf"$\mathregular{{{int(x):,}}}$"))

# ax.tick_params(axis="y", labelcolor=C_temp, labelsize=9, which="both")
# ax.tick_params(axis="x", labelsize=9)
# plt.setp(ax.get_yticklabels(which="both"), color=C_temp)

# # ── Density ───────────────────────────────────────────────────────────────────
# ax2 = ax.twinx()
# ax2.spines["top"].set_visible(True)

# ax2.plot(radius, rho_cgs,
#          color=C_rho, lw=2.2, ls=(0,(6,2)),
#          solid_capstyle="round", zorder=5)

# ax2.set_ylabel(r"Density  [g cm$^{-3}$]", fontsize=11, color=C_rho, labelpad=6)
# ax2.tick_params(axis="y", labelcolor=C_rho, direction="out", labelsize=9)

# ax2.set_yscale("log")
# ax2.set_ylim([rho_cgs.min()*0.4, rho_cgs.max()*4])

# ax2.yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))

# # ── Legend ────────────────────────────────────────────────────────────────────
# leg_handles = [
#     mpl.lines.Line2D([0],[0], color=C_temp, lw=2.2,           label="Temperature"),
#     mpl.lines.Line2D([0],[0], color=C_rho,  lw=2.2, ls=(0,(6,2)), label="Density"),
# ]

# ax.legend(handles=leg_handles,
#           loc="lower center", bbox_to_anchor=(0.5, 0.02),
#           fontsize=9.5, framealpha=0.95,
#           edgecolor="#aaaaaa", handlelength=2.8,
#           borderpad=0.7)

# # plt.savefig("modelS_temp_density.pdf", bbox_inches="tight")
# # plt.savefig("solar_model_s_temp_density_interior.pdf", bbox_inches="tight")

# plt.show()
# # %%
# # %%
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.ticker import LogFormatterSciNotation, AutoMinorLocator
# from scipy.constants import k, m_u

# # --- Load VAL3C atmosphere ---
# height, T_atm, pe_arr, pg = np.loadtxt(
#     "VAL3C_rutten.dat",
#     unpack=True,
#     skiprows=1,
#     usecols=(0,3,4,5)
# )

# pg /= 10
# pe_arr /= 10

# mu = 1.0/(1.0 + pe_arr/pg)
# rho = (mu*m_u*pg)/(k*T_atm)
# rho_atm = rho*1e-3

# # --- Load solar interior ---
# model = np.genfromtxt(
#     "https://users-phys.au.dk/~jcd/solar_models/cptrho.l5bi.d.15c",
#     comments="#"
# )

# r = model[:,0]
# rho_int = model[:,2]
# T_int = model[:,5]

# # ── Style (UNCHANGED from your code) ─────────────────────────────────────────
# plt.rcParams.update({
#     "font.family": "serif",
#     "mathtext.fontset": "cm",
#     "axes.linewidth": 1.1,
#     "xtick.direction": "out",
#     "ytick.direction": "out",
#     "xtick.major.size": 6,
#     "ytick.major.size": 6,
#     "xtick.minor.size": 3.5,
#     "ytick.minor.size": 3.5,
# })

# # Atmosphere region colours
# C_photo  = "#2a7d4f"
# C_chrom  = "#b5480f"
# C_trans  = "#7b2d8b"
# C_corona = "#a6122d"

# # Interior region colours
# C_core = "#e6b800"
# C_rad  = "#c85c00"
# C_cz   = "#6c8fbf"

# # Data colours (same as your plot)
# C_temp = "#1464a0"
# C_rho  = "#c47a00"

# # --- Figure ---
# fig,(ax1,ax2)=plt.subplots(
#     1,2,
#     figsize=(15,5)
# )

# fig.subplots_adjust(
#     left=0.08,
#     right=0.88,
#     top=0.93,
#     bottom=0.13,
#     wspace=0.30
# )

# # =========================================================
# # (a) SOLAR INTERIOR
# # =========================================================

# regions_int=[
#     (0,0.25,C_core,"Core"),
#     (0.25,0.71,C_rad,"Radiative\nZone"),
#     (0.71,1.0,C_cz,"Convection\nZone"),
# ]

# for x0,x1,color,label in regions_int:

#     ax1.axvspan(x0,x1,alpha=0.18,color=color,lw=0)

#     if x0!=0:
#         ax1.axvline(x=x0,color=color,lw=1.0,ls=(0,(5,3)),alpha=0.7)

#     ax1.text(
#         (x0+x1)/2,
#         4e6,
#         label,
#         ha="center",
#         va="center",
#         fontsize=12,
#         color=color,
#         fontweight="bold"
#     )

# ax1.plot(r,T_int,color=C_temp,lw=2.2)

# ax1.set_xlabel(r"Radius [$R_\odot$]", fontsize=11)
# ax1.set_ylabel("Temperature [K]",color=C_temp, fontsize=11)
# ax1.set_yscale("log")

# ax1.tick_params(axis="y",labelcolor=C_temp,labelsize=9)
# ax1.tick_params(axis="x",labelsize=9)

# ax1b=ax1.twinx()

# ax1b.plot(
#     r,
#     rho_int,
#     color=C_rho,
#     lw=2.2,
#     ls=(0,(6,2))
# )

# ax1b.set_ylabel(r"Density [g cm$^{-3}$]",color=C_rho, fontsize=11)
# ax1b.set_yscale("log")
# ax1b.tick_params(axis="y",labelcolor=C_rho,labelsize=9)

# ax1b.yaxis.set_major_formatter(
#     LogFormatterSciNotation(labelOnlyBase=False)
# )

# ax1.text(
#     -0.12,1.03,"(a)",
#     transform=ax1.transAxes,
#     fontsize=16,
#     fontweight="bold"
# )

# # =========================================================
# # (b) ATMOSPHERE
# # =========================================================

# regions_atm=[
#     (height.min(),500,C_photo,"Photosphere",False),
#     (500,2100,C_chrom,"Chromosphere",False),
#     (2100,2300,C_trans,"Transition\nRegion",True),
#     (2300,height.max(),C_corona,"Corona",True),
# ]

# for x0,x1,color,label,vert in regions_atm:

#     ax2.axvspan(x0,x1,alpha=0.18,color=color,lw=0)

#     if x0!=height.min():
#         ax2.axvline(x=x0,color=color,lw=1.0,ls=(0,(5,3)),alpha=0.7)

#     ax2.text(
#         (x0+x1)/2,
#         25000,
#         label,
#         ha="center",
#         va="center",
#         fontsize=12,
#         color=color,
#         rotation=90 if vert else 0,
#         fontweight="bold"
#     )

# ax2.plot(height,T_atm,color=C_temp,lw=2.2)

# ax2.set_xlabel("Height [km]", fontsize=11)
# ax2.set_ylabel("Temperature [K]",color=C_temp, fontsize=11)

# ax2.set_yscale("log")
# ax2.set_ylim(ax2.get_ylim()[0], 1e5)

# ax2.tick_params(axis="y",labelcolor=C_temp,labelsize=9)
# ax2.tick_params(axis="x",labelsize=9)

# ax2b=ax2.twinx()

# ax2b.plot(
#     height,
#     rho_atm,
#     color=C_rho,
#     lw=2.2,
#     ls=(0,(6,2))
# )

# ax2b.set_ylabel(r"Density [g cm$^{-3}$]",color=C_rho, fontsize=11)
# ax2b.set_yscale("log")

# ax2b.tick_params(axis="y",labelcolor=C_rho,labelsize=9)

# ax2b.yaxis.set_major_formatter(
#     LogFormatterSciNotation(labelOnlyBase=False)
# )

# ax2.text(
#     -0.12,1.03,"(b)",
#     transform=ax2.transAxes,
#     fontsize=12,
#     fontweight="bold"
# )

# # =========================================================
# # Shared legend
# # =========================================================

# handles=[
#     mpl.lines.Line2D([0],[0],color=C_temp,lw=2.2,label="Temperature"),
#     mpl.lines.Line2D([0],[0],color=C_rho,lw=2.2,ls=(0,(6,2)),label="Density"),
# ]

# fig.legend(
#     handles=handles,
#     loc="lower center",
#     ncol=2,
#     bbox_to_anchor=(0.5,0.00),
#     fontsize=12,
#     framealpha=0.95,
#     edgecolor="#aaaaaa"
# )

# plt.savefig("solar_structure_two_panels.pdf",bbox_inches="tight")
# plt.show()
# # %%


# %%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation
from scipy.constants import k, m_u
from matplotlib.gridspec import GridSpec

# --- Load VAL3C atmosphere ---
height, T_atm, pe_arr, pg = np.loadtxt(
    "VAL3C_rutten.dat",
    unpack=True,
    skiprows=1,
    usecols=(0,3,4,5)
)

pg /= 10
pe_arr /= 10

mu = 1.0/(1.0 + pe_arr/pg)
rho = (mu*m_u*pg)/(k*T_atm)
rho_atm = rho*1e-3

# --- Load solar interior ---
model = np.genfromtxt(
    "https://users-phys.au.dk/~jcd/solar_models/cptrho.l5bi.d.15c",
    comments="#"
)

r = model[:,0]
rho_int = model[:,2]
T_int = model[:,5]

# ── Style (bigger fonts) ─────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.linewidth": 1.2,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "font.size": 13,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# Colors
C_photo  = "#2a7d4f"
C_chrom  = "#b5480f"
C_trans  = "#7b2d8b"
C_corona = "#a6122d"

C_core = "#e6b800"
C_rad  = "#c85c00"
C_cz   = "#6c8fbf"
C_tacho = "#444444"  # tachocline (neutral dark)

C_temp = "#1464a0"
C_rho  = "#c47a00"

# =========================================================
# FIGURE (VERTICAL LAYOUT)
# =========================================================

fig = plt.figure(figsize=(8,10))
gs = GridSpec(2, 1, height_ratios=[1,1], hspace=0.25)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# =========================================================
# (a) SOLAR INTERIOR
# =========================================================

regions_int = [
    (0,0.2,C_core,"Core"),
    (0.2,0.71,C_rad,"Radiative\nZone"),
    (0.71,1.0,C_cz,"Convection\nZone"),
]

for x0,x1,color,label in regions_int:
    ax1.axvspan(x0,x1,alpha=0.18,color=color,lw=0)

    if x0!=0:
        ax1.axvline(x=x0,color=color,lw=1.2,ls=(0,(5,3)),alpha=0.7)

    ax1.text((x0+x1)/2,4e6,label,
             ha="center",va="center",
             fontsize=13,color=color,fontweight="bold")

# --- Tachocline (NEW) ---
tacho_min, tacho_max = 0.69, 0.72
ax1.axvspan(tacho_min, tacho_max, color=C_tacho, alpha=0.25, zorder=1)
ax1.text(0.705, 1e6, "Tachocline",
         ha="center", va="center",
         fontsize=12, color=C_tacho, rotation=90)

# Data
ax1.plot(r,T_int,color=C_temp,lw=2.5)

ax1.set_xlabel(r"Radius [$R_\odot$]")
ax1.set_ylabel("Temperature [K]",color=C_temp)
ax1.set_yscale("log")

ax1.tick_params(axis="y",labelcolor=C_temp)

ax1b = ax1.twinx()
ax1b.plot(r,rho_int,color=C_rho,lw=2.5,ls=(0,(6,2)))

ax1b.set_ylabel(r"Density [g cm$^{-3}$]",color=C_rho)
ax1b.set_yscale("log")
ax1b.tick_params(axis="y",labelcolor=C_rho)

ax1b.yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))

ax1.text(-0.1,1.03,"(a)",transform=ax1.transAxes,
         fontsize=16,fontweight="bold")

# =========================================================
# (b) ATMOSPHERE
# =========================================================

regions_atm = [
    (height.min(),500,C_photo,"Photosphere",False),  # UPDATED
    (500,2100,C_chrom,"Chromosphere",False),
    (2100,2300,C_trans,"Transition\nRegion",True),
    (2300,height.max(),C_corona,"Corona",True),
]

for x0,x1,color,label,vert in regions_atm:
    ax2.axvspan(x0,x1,alpha=0.18,color=color,lw=0)

    if x0!=height.min():
        ax2.axvline(x=x0,color=color,lw=1.2,ls=(0,(5,3)),alpha=0.7)

    ax2.text((x0+x1)/2,25000,label,
             ha="center",va="center",
             fontsize=13,color=color,
             rotation=90 if vert else 0,
             fontweight="bold")

# Data
ax2.plot(height,T_atm,color=C_temp,lw=2.5)

ax2.set_xlabel("Height [km]")
ax2.set_ylabel("Temperature [K]",color=C_temp)

ax2.set_yscale("log")
ax2.set_ylim(ax2.get_ylim()[0],1e5)

ax2.tick_params(axis="y",labelcolor=C_temp)

ax2b = ax2.twinx()
ax2b.plot(height,rho_atm,color=C_rho,lw=2.5,ls=(0,(6,2)))

ax2b.set_ylabel(r"Density [g cm$^{-3}$]",color=C_rho)
ax2b.set_yscale("log")
ax2b.tick_params(axis="y",labelcolor=C_rho)

ax2b.yaxis.set_major_formatter(LogFormatterSciNotation(labelOnlyBase=False))

ax2.text(-0.1,1.03,"(b)",transform=ax2.transAxes,
         fontsize=16,fontweight="bold")

# =========================================================
# LEGEND
# =========================================================

handles = [
    mpl.lines.Line2D([0],[0],color=C_temp,lw=2.5,label="Temperature"),
    mpl.lines.Line2D([0],[0],color=C_rho,lw=2.5,ls=(0,(6,2)),label="Density"),
]

fig.legend(handles=handles,
           loc="lower center",
           ncol=2,
           bbox_to_anchor=(0.5,0.00),
           fontsize=13,
           framealpha=0.95)

plt.savefig("solar_structure_vertical.pdf",bbox_inches="tight")
plt.show()
# %%
