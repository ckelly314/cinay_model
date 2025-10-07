import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Figure Settings
dpi = 500
img_format = "png"
fs = 15
ff = "arial"


# Load Falkor Bottle Data
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="sigma0")
st9 = falkor.loc[falkor["Station"] == 9]
st9 = st9.loc[(st9["sigma0"] > 25.0) & (st9["sigma0"] < 27.4)]
st2 = falkor.loc[falkor["Station"] == 2]
st2 = st2.loc[(st2["sigma0"] > 25.0) & (st2["sigma0"] < 27.4)]


# Plot Tracer Profiles
fig, ax = plt.subplots(1, 7, figsize=(12, 5))
tracers = list(("O2", "NO3", "NO2", "NH4", "DIP", "TA", "pH"))
tracer_labels = list(
    (
        "$[O_2]$ ($\mu$$mol$ $kg^{-1})$",
        "$[NO_3^-]$ ($\mu$$mol$ $kg^{-1})$",
        "$[NO_2^-]$ ($\mu$$mol$ $kg^{-1})$",
        "$[NH_4^+]$ ($\mu$$mol$ $kg^{-1})$",
        "$[DIP]$ ($\mu$$mol$ $kg^{-1})$",
        "$ TA$ ($\mu$$mol$ $kg^{-1})$",
        "$pH$",
    )
)
lowlims = np.array((0, 20, -0.03, -0.001, 2, 2270, 7.54))
uplims = np.array((15, 45, 3, 0.06, 3.6, 2380, 7.62))
colors = list(
    (
        "mediumvioletred",
        "blueviolet",
        "royalblue",
        "forestgreen",
        "gold",
        "sandybrown",
        "tomato",
    )
)
for i in np.arange(0, 7):
    ax[i].plot(
        st9[tracers[i]],
        st9["sigma0"],
        linewidth=2,
        ls="-",
        marker="o",
        markersize=5,
        fillstyle="none",
        color=colors[i],
    )
    ax[i].set_ylim((25, 27.4))
    ax[i].invert_yaxis()
    ax[i].xaxis.set_label_position("top")
    ax[i].xaxis.set_ticks_position("top")
    ax[i].tick_params(axis="x", labelsize=fs)

    if i == 0:
        ax[i].set_ylabel(
            "$\sigma$$_{\\theta}$ $(kg$ $m^{-3}$)", fontfamily=ff, fontsize=fs - 2
        )
        ax[i].tick_params(axis="y", labelsize=fs)
    else:
        ax[i].axes.yaxis.set_ticklabels([])

    ax[i].set_xlim((lowlims[i], uplims[i]))
    ax[i].set_xlabel(tracer_labels[i], fontfamily=ff, fontsize=fs - 3)


fig.tight_layout()
plt.savefig("figures/station_plots/station9.{}".format(img_format), dpi=dpi)
