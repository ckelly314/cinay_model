## Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sc
import statsmodels.api as sm

# Set Directory
file_path = "OM_variations/anderson/{}"
data_path = "output/{}".format(file_path)
fig_path = "figures/{}"  # .format(file_path)
fig_format = "pdf"
dpi = 500
fs = 13
ff = "arial"

# stations = np.arange(1,15,1, dtype=float)
stations = np.array((1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19))

# Import Parameters and Results
falkor = pd.read_csv("falkor_clean.csv")
coeffs_mean = np.array(pd.read_csv(data_path.format("coeffs_mean.csv")))
coeffs_se = np.array(pd.read_csv(data_path.format("coeffs_se.csv")))
slopes_mean = np.array(pd.read_csv(data_path.format("slopes_mean.csv")))
slopes_se = np.array(pd.read_csv(data_path.format("slopes_se.csv")))
residuals_perc = np.array(pd.read_csv(data_path.format("percentage_residuals.csv")))
residuals = np.array(pd.read_csv(data_path.format("residuals.csv")))
relimp_mean = np.array(pd.read_csv(data_path.format("relative_importances_mean.csv")))
relimp_se = np.array(pd.read_csv(data_path.format("relative_importances_se.csv"))) * 1.0
relimp_nit = np.array(
    pd.read_csv(data_path.format("relative_importances_nitrite_mean.csv"))
)
relimp_nit_se = (
    np.array(pd.read_csv(data_path.format("relative_importances_nitrite_se.csv"))) * 1.0
)
layers = np.array(pd.read_csv(data_path.format("layers.csv")))

# Layer Labels
ylab = list()
for i in np.arange(0, len(layers) - 1):
    if i == 0:
        ylab.append("Layer {}".format(int(i + 1)))
    else:
        ylab.append("{}".format(int(i + 1)))

# Coefficients Heat Maps
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
hm = sns.heatmap(
    ax=ax1,
    data=coeffs_mean,
    annot=True,
    annot_kws={"size": 7},
    fmt=".4f",
    cmap="Spectral_r",
    vmin=0,
    vmax=1,
    linewidths=0.5,
)
hm.set_yticks(ticks=np.arange(0.5, 16.5, 1))
hm.set_yticklabels(labels=ylab, size=8, rotation=0)
hm.set_xticklabels(
    ["DNRN", "Denitrification", "Anammox", "Nitrite Oxidation", "$CaCO_3$ Dissolution"],
    rotation=45,
    horizontalalignment="right",
)
hm.set_title("Mean Relative Reaction Rates")

hm = sns.heatmap(
    ax=ax2,
    data=coeffs_se,
    annot=True,
    annot_kws={"size": 7},
    fmt=".4f",
    cmap="Spectral_r",
    vmin=0,
    vmax=0.1,
    linewidths=0.5,
)
hm.set_yticks(ticks=np.arange(0.5, 16.5, 1))
hm.set_yticklabels(labels=ylab, size=8, rotation=0)
hm.set_xticklabels(
    ["DNRN", "Denitrification", "Anammox", "Nitrite Oxidation", "$CaCO_3$ Dissolution"],
    rotation=45,
    horizontalalignment="right",
)
hm.set_title("Standard Deviation")
fig1.tight_layout(pad=2.5)
plt.savefig(
    fig_path.format("reaction_coefficients.{}".format(fig_format))
)  # , dpi = dpi)

# Residuals  Heat Maps
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
hm = sns.heatmap(
    ax=ax1,
    data=residuals,
    annot=True,
    annot_kws={"size": 6},
    fmt=".4f",
    cmap="Spectral_r",
    vmin=-0.10,
    vmax=0.10,
    linewidths=0.5,
)
hm.set_yticks(ticks=np.arange(0.5, 16.5, 1))
hm.set_yticklabels(labels=ylab, size=8, rotation=0)
hm.set_xticklabels(
    [
        "$\Delta$$NO_3^-$",
        "$\Delta$$NO_2^-$",
        "$\Delta$$NH_4^+$",
        "$\Delta$$N^*$",
        "$\Delta$$TA$",
        "$\Delta$$DIC$",
    ]
)
hm.set_title("Residuals")


hm = sns.heatmap(
    ax=ax2,
    data=residuals_perc,
    annot=True,
    annot_kws={"size": 6},
    cmap="Spectral_r",
    vmin=0,
    vmax=100,
    linewidths=0.5,
)
hm.set_yticks(ticks=np.arange(0.5, 16.5, 1))
hm.set_yticklabels(labels=ylab, size=8, rotation=0)
hm.set_xticklabels(
    [
        "$\Delta$$NO_3^-$",
        "$\Delta$$NO_2^-$",
        "$\Delta$$NH_4^+$",
        "$\Delta$$N^*$",
        "$\Delta$$TA$",
        "$\Delta$$DIC$",
    ]
)
hm.set_title("Percentage Residuals")
fig2.tight_layout(pad=2.5)
plt.savefig(fig_path.format("residuals.{}".format(fig_format)), dpi=dpi)

# Relative Contributions
new2_dnrn = relimp_nit[:, 3] / (relimp_nit[:, 1] + relimp_nit[:, 3]) * 100
new2_denit = (relimp_nit[:, 1]) / (relimp_nit[:, 1] + relimp_nit[:, 3]) * 100
new2_se = (
    (relimp_nit_se[:, 3] / relimp_nit[:, 3])
    + (relimp_nit_se[:, 1] + relimp_nit_se[:, 3])
    / (relimp_nit[:, 1] + relimp_nit[:, 3])
) * new2_dnrn
new3_prod = (
    relimp_nit[:, 3]
    / (relimp_nit[:, 0] + relimp_nit[:, 1] + relimp_nit[:, 2] + relimp_nit[:, 3])
    * 100
)
new3_cons = (
    (relimp_nit[:, 0] + relimp_nit[:, 1] + relimp_nit[:, 2])
    / (relimp_nit[:, 0] + relimp_nit[:, 1] + relimp_nit[:, 2] + relimp_nit[:, 3])
    * 100
)
new3_se = (
    (relimp_nit_se[:, 3] / relimp_nit[:, 3])
    + (
        relimp_nit_se[:, 0]
        + relimp_nit_se[:, 1]
        + relimp_nit_se[:, 2]
        + relimp_nit_se[:, 3]
    )
    / (relimp_nit[:, 0] + relimp_nit[:, 1] + relimp_nit[:, 2] + relimp_nit[:, 3])
) * new3_prod

fig3, ((ax1, ax2, ax5), (ax3, ax6, ax4)) = plt.subplots(
    nrows=2, ncols=3, figsize=(15, 10)
)
ax1.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 0],
    xerr=relimp_se[:, 0],
    label="Anammox",
    color="thistle",
)
ax1.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 1],
    left=relimp_mean[:, 0],
    label="Denitrification",
    color="wheat",
)
xx = np.linspace(0, 110, 10)
yy = np.ones(np.shape(xx))
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = relimp_mean[i, 0] * np.ones(np.shape(yy1))
    ax1.plot(xx1, yy1, color="black", linewidth=1.0)
ax1.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax1.set(xlim=(0, 105))
ax1.set_yticks(ticks=np.arange(1, 17))
ax1.set_yticklabels(labels=ylab, size=16)
ax1.tick_params(axis="x", labelsize=fs + 3)
ax1.tick_params(axis="y", labelsize=fs + 3)
ax1.axes.xaxis.set_ticklabels([])
ax1.legend(loc=3, fontsize=13)
ax1.invert_yaxis()

ax2.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 2],
    xerr=relimp_se[:, 2],
    label="Nitrite Oxidation",
    color="thistle",
)
ax2.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 3],
    left=relimp_mean[:, 2],
    label="Nitrite Reduction",
    color="wheat",
)
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = relimp_mean[i, 2] * np.ones(np.shape(yy1))
    ax2.plot(xx1, yy1, color="black", linewidth=1.0)
ax2.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax2.set(xlim=(0, 105))
ax2.set_yticks(ticks=np.arange(1, 17))
ax2.set_yticklabels(labels=ylab, size=16)
ax2.tick_params(axis="x", labelsize=fs + 3)
ax2.tick_params(axis="y", labelsize=fs + 3)
ax2.axes.xaxis.set_ticklabels([])
ax2.axes.yaxis.set_ticklabels([])
ax2.legend(loc=3, fontsize=13)
ax2.invert_yaxis()

ax3.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 5],
    xerr=relimp_se[:, 5],
    label="DNRN",
    color="thistle",
)
ax3.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 4],
    left=relimp_mean[:, 5],
    label="Nitrite Oxidation",
    color="wheat",
)
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = relimp_mean[i, 5] * np.ones(np.shape(yy1))
    ax3.plot(xx1, yy1, color="black", linewidth=1.0)
ax3.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax3.set(xlim=(0, 105))
ax3.set_xlabel("Relative Contribution %", fontsize=fs + 3)
ax3.set_yticks(ticks=np.arange(1, 17))
ax3.set_yticklabels(labels=ylab, size=16)
ax3.tick_params(axis="x", labelsize=fs + 3)
ax3.tick_params(axis="y", labelsize=fs + 3)
ax3.legend(loc=3, fontsize=13)
ax3.invert_yaxis()

ax4.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 7],
    xerr=relimp_se[:, 7],
    label="Other Reactions",
    color="thistle",
)
ax4.barh(
    y=np.arange(1, len(layers)),
    width=relimp_mean[:, 6],
    left=relimp_mean[:, 7],
    label="$CaCO_3$ Dissolution",
    color="wheat",
)
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = relimp_mean[i, 7] * np.ones(np.shape(yy1))
    ax4.plot(xx1, yy1, color="black", linewidth=1.0)
ax4.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax4.set(xlim=(0, 105))
ax4.set_xlabel("Relative Contribution %", fontsize=fs + 3)
ax4.set_yticks(ticks=np.arange(1, 17))
ax4.set_yticklabels(labels=ylab, size=16)
ax4.tick_params(axis="x", labelsize=fs + 3)
ax4.tick_params(axis="y", labelsize=fs + 3)
ax4.axes.yaxis.set_ticklabels([])
ax4.legend(loc=3, fontsize=13)
ax4.invert_yaxis()

ax5.barh(
    y=np.arange(1, len(layers)),
    width=new2_dnrn,
    xerr=new2_se,
    label="DNRN",
    color="thistle",
)
ax5.barh(
    y=np.arange(1, len(layers)),
    width=new2_denit,
    left=new2_dnrn,
    label="Denitrification",
    color="wheat",
)
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = new2_dnrn[i] * np.ones(np.shape(yy1))
    ax5.plot(xx1, yy1, color="black", linewidth=1.0)
ax5.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax5.set(xlim=(0, 105))
ax5.set_yticks(ticks=np.arange(1, 17))
ax5.set_yticklabels(labels=ylab, size=16)
ax5.legend(loc=3, fontsize=13)
ax5.tick_params(axis="x", labelsize=fs + 3)
ax5.tick_params(axis="y", labelsize=fs + 3)
ax5.axes.xaxis.set_ticklabels([])
ax5.axes.yaxis.set_ticklabels([])
ax5.invert_yaxis()

ax6.barh(
    y=np.arange(1, len(layers)),
    width=new3_prod,
    xerr=new3_se,
    label="Nitrite Production",
    color="thistle",
)
ax6.barh(
    y=np.arange(1, len(layers)),
    width=new3_cons,
    left=new3_prod,
    label="Nitrite Consumption",
    color="wheat",
)
for i in np.arange(0, 16):
    yy1 = np.linspace((i + 1) - 0.36, (i + 1) + 0.36, 5)
    xx1 = new3_prod[i] * np.ones(np.shape(yy1))
    ax6.plot(xx1, yy1, color="black", linewidth=1.0)
ax6.axvline(x=50, color="red", linewidth=1.5, ls="--")
ax6.set(xlim=(0, 105))
ax6.set_xlabel("Relative Contribution %", fontsize=fs + 3)
ax6.set_yticks(ticks=np.arange(1, 17))
ax6.invert_yaxis()
ax6.set_yticklabels(labels=ylab, size=16)
ax6.tick_params(axis="x", labelsize=fs + 3)
ax6.tick_params(axis="y", labelsize=fs + 3)
ax6.axes.yaxis.set_ticklabels([])
ax6.legend(loc=3, fontsize=13)


fig3.tight_layout(pad=2.5)
plt.savefig(fig_path.format("relative_importances.{}".format(fig_format)), dpi=dpi)
