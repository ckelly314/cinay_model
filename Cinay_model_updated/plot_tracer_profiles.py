## Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sc

# Figure Settings
fig_path = "figures/station_plots/"
ftype = ".png"
dpi = 500
fs = 20
ff = "arial"

# Import Parameters and Results
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="sigma0")
layers = np.array(pd.read_csv("output/initial_results/layers.csv"))
stations = np.array((1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19))
falkor_CTD = pd.read_csv(
    "../../../ETP Storms/Data/FK180624/Falkor_CTD_casts/FK180624_CTD.csv"
)

# Select Correct Station
falkor = falkor.iloc[np.where(np.isin(falkor["Station"], stations))].sort_values(
    by=["Station", "sigma0"]
)
sta = falkor.loc[falkor["Station"] == 2]
stb = falkor.loc[falkor["Station"] == 7]
st9 = falkor.loc[falkor["Station"] == 9]
stc = falkor.loc[falkor["Station"] == 12]
std = falkor.loc[falkor["Station"] == 16]

# Select CTD O2
st9_ctd = falkor_CTD.loc[(falkor_CTD["station"] == 9.0) & (falkor_CTD["cast"] == 1.0)]
sta_ctd = falkor_CTD.loc[(falkor_CTD["station"] == 2.0) & (falkor_CTD["cast"] == 1.0)]
stb_ctd = falkor_CTD.loc[(falkor_CTD["station"] == 7.0) & (falkor_CTD["cast"] == 1.0)]
stc_ctd = falkor_CTD.loc[(falkor_CTD["station"] == 12.0) & (falkor_CTD["cast"] == 2.0)]
std_ctd = falkor_CTD.loc[(falkor_CTD["station"] == 16.0) & (falkor_CTD["cast"] == 1.0)]

# Depth Profiles (use CTD data)
fig3, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 7))
sns.despine(fig=fig3, right=True, top=False, bottom=True)
# ax1.scatter(falkor['O2'], falkor['sigma0'], s = 3)
ax1.plot(
    sta["O2"],
    sta["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumturquoise",
    label="Station 2",
    alpha=0.5,
    zorder=1,
)
ax1.plot(
    stb["O2"],
    stb["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumorchid",
    label="Station 7",
    alpha=0.5,
    zorder=2,
)
ax1.plot(
    st9["O2"],
    st9["sigma0"],
    linewidth=2.0,
    marker="o",
    markersize=5,
    fillstyle="none",
    color="black",
    label="Station 9",
    zorder=5,
)
ax1.plot(
    stc["O2"],
    stc["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="lightcoral",
    label="Station 12",
    alpha=0.5,
    zorder=3,
)
ax1.plot(
    std["O2"],
    std["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="springgreen",
    label="Station 16",
    alpha=0.5,
    zorder=4,
)

ax1.set_ylim((25, 27.4))
ax1.set_xlim((0, 15))
ax1.invert_yaxis()
ax1.set_xlabel("$[O_2]$ ($\mu$$mol$ $kg^{-1}) $", fontfamily=ff, fontsize=fs)
ax1.set_ylabel("$\sigma$$_{\\theta}$ ($kg$ $m^{-3}$)", fontfamily=ff, fontsize=fs)
ax1.tick_params(axis="x", labelsize=fs)
ax1.tick_params(axis="y", labelsize=fs)

ax1.xaxis.set_label_position("top")
ax1.xaxis.set_ticks_position("top")
ax1.xaxis.set_label_coords(0.5, 1.1)
ax1.legend(loc=4, fontsize=fs - 5)
for i in np.arange(0, len(layers) - 1):
    x = np.arange(0, 16, 1)
    if (i % 2) == 0:
        ax1.fill_between(x, layers[i], layers[i + 1], color="whitesmoke", alpha=0.4)
    elif (i % 2) != 0:
        ax1.fill_between(x, layers[i], layers[i + 1], color="gainsboro", alpha=0.4)

# ax2.scatter(falkor['NO2'], falkor['sigma0'], s= 3)
ax2.plot(
    st9["NO2"],
    st9["sigma0"],
    linewidth=2.0,
    marker="o",
    markersize=5,
    fillstyle="none",
    color="black",
    label="Station 9",
    zorder=5,
)
ax2.plot(
    sta["NO2"],
    sta["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumturquoise",
    label="Station 2",
    alpha=0.5,
    zorder=1,
)
ax2.plot(
    stb["NO2"],
    stb["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumorchid",
    label="Station 7",
    alpha=0.5,
    zorder=2,
)
ax2.plot(
    stc["NO2"],
    stc["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="lightcoral",
    label="Station 12",
    alpha=0.5,
    zorder=3,
)
ax2.plot(
    std["NO2"],
    std["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="springgreen",
    label="Station 16",
    alpha=0.5,
    zorder=4,
)

ax2.set_ylim((25, 27.4))
ax2.set_xlim((0, 3))
ax2.invert_yaxis()
ax2.set_xlabel("$[NO_2^-]$ ($\mu$$mol$ $kg^{-1}) $", fontfamily=ff, fontsize=fs)
ax2.axes.yaxis.set_ticklabels([])
# ax2.set_ylabel('$\sigma$$_0$ (kg/$m^3$)')
ax2.tick_params(axis="x", labelsize=fs)
ax2.tick_params(axis="y", labelsize=fs)
ax2.xaxis.set_label_position("top")
ax2.xaxis.set_ticks_position("top")
ax2.xaxis.set_label_coords(0.5, 1.1)

for i in np.arange(0, len(layers) - 1):
    # if i == 0:
    #    txt = 'Layer 1'
    # else:
    #    txt = '{}'.format(i+1)
    # ax2.text(2.5, (layers[i]+ layers[i+1])/2+0.01, txt, fontsize =7 )
    x = np.arange(0, 4, 1)
    if (i % 2) == 0:
        ax2.fill_between(x, layers[i], layers[i + 1], color="whitesmoke", alpha=0.4)
    elif (i % 2) != 0:
        ax2.fill_between(x, layers[i], layers[i + 1], color="gainsboro", alpha=0.4)


# ax3.scatter(falkor['pH'], falkor['sigma0'], s=3)ax1.plot(st9['O2'], st9['sigma0'], linewidth = 2.0, marker = 'o', markersize = 5, fillstyle = 'none', color = 'black', label = 'Station 9', zorder = 5)
ax3.plot(
    st9["pH"],
    st9["sigma0"],
    linewidth=2.0,
    marker="o",
    markersize=5,
    fillstyle="none",
    color="black",
    label="Station 9",
    zorder=5,
)
ax3.plot(
    sta["pH"],
    sta["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumturquoise",
    label="Station 2",
    alpha=0.5,
    zorder=1,
)
ax3.plot(
    stb["pH"],
    stb["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="mediumorchid",
    label="Station 7",
    alpha=0.5,
    zorder=2,
)
ax3.plot(
    stc["pH"],
    stc["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="lightcoral",
    label="Station 12",
    alpha=0.5,
    zorder=3,
)
ax3.plot(
    std["pH"],
    std["sigma0"],
    linewidth=1.0,
    marker="o",
    markersize=2,
    fillstyle="none",
    color="springgreen",
    label="Station 16",
    alpha=0.5,
    zorder=4,
)

ax3.set_ylim((25, 27.4))
ax3.set_xlim((7.54, 7.65))
ax3.invert_yaxis()
ax3.set_xlabel("pH", fontfamily=ff, fontsize=fs)
ax3.axes.yaxis.set_ticklabels([])
ax3.tick_params(axis="x", labelsize=fs)
ax3.tick_params(axis="y", labelsize=fs)
# ax3.set_ylabel('$\sigma$$_0$ (kg/$m^3$)')
ax3.xaxis.set_label_position("top")
ax3.xaxis.set_ticks_position("top")
ax3.xaxis.set_label_coords(0.5, 1.1)
for i in np.arange(0, len(layers) - 1):
    x = np.arange(7.54, 7.66, 0.01)
    if (i % 2) == 0:
        ax3.fill_between(x, layers[i], layers[i + 1], color="whitesmoke", alpha=0.4)
    elif (i % 2) != 0:
        ax3.fill_between(x, layers[i], layers[i + 1], color="gainsboro", alpha=0.4)

fig3.tight_layout(pad=1.5)
plt.savefig("{}tracer_depth_profiles{}".format(fig_path, ftype), dpi=dpi)
