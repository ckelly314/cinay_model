## Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set Path
data_path = "output/OM_variations/{}/{}"
fig_path = "figures/OM_variations/{}/{}"
fig_format = "png"
dpi = 500
fs = 13
ff = "arial"

# Folder Names
Cox_list = np.arange(-4.0, 4.5, 0.5)
fpath = list()
for i in Cox_list:
    if i < 0:
        fpath.append("cox_neg{}".format(round(abs(i), 2)))
    elif i == 0:
        fpath.append("cox_zero")
    elif i > 0:
        fpath.append("cox_pos{}".format(round(i, 1)))

# Get Relative Importance Data
RI_anmx = np.zeros((len(Cox_list), 6))
RI_denit = np.zeros((len(Cox_list), 6))
RI_nitox = np.zeros((len(Cox_list), 6))
RI_DNRN = np.zeros((len(Cox_list), 6))
for i in np.arange(0, len(Cox_list)):
    nit_sys = np.array(
        pd.read_csv(data_path.format(fpath[i], "relative_importances_nitrite_mean.csv"))
    )
    RI_anmx[i, :] = np.array(
        [
            np.mean(nit_sys[:, 0]),
            nit_sys[3, 0],
            nit_sys[5, 0],
            nit_sys[7, 0],
            nit_sys[10, 0],
            nit_sys[13, 0],
        ]
    )
    RI_denit[i, :] = np.array(
        [
            np.mean(nit_sys[:, 1]),
            nit_sys[3, 1],
            nit_sys[5, 1],
            nit_sys[7, 1],
            nit_sys[10, 1],
            nit_sys[13, 1],
        ]
    )
    RI_nitox[i, :] = np.array(
        [
            np.mean(nit_sys[:, 2]),
            nit_sys[3, 2],
            nit_sys[5, 2],
            nit_sys[7, 2],
            nit_sys[10, 2],
            nit_sys[13, 2],
        ]
    )
    RI_DNRN[i, :] = np.array(
        [
            np.mean(nit_sys[:, 3]),
            nit_sys[3, 3],
            nit_sys[5, 3],
            nit_sys[7, 3],
            nit_sys[10, 3],
            nit_sys[13, 3],
        ]
    )
    relimps = np.array(
        pd.read_csv(data_path.format(fpath[i], "relative_importances_mean.csv"))
    )


# Plot Ratios
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 3))
# 'Denitrification/Anammox'
ax1.plot(
    Cox_list,
    RI_anmx[:, 1] / (RI_anmx[:, 1] + RI_denit[:, 1]) * 100,
    color="mediumturquoise",
    ls="solid",
    label="ODZ Top (Layer 4)",
    linewidth=1,
)
ax1.plot(
    Cox_list,
    RI_anmx[:, 2] / (RI_anmx[:, 2] + RI_denit[:, 2]) * 100,
    color="springgreen",
    ls="solid",
    label="SNM (Layer 6)",
    linewidth=1,
)
ax1.plot(
    Cox_list,
    RI_anmx[:, 3] / (RI_anmx[:, 3] + RI_denit[:, 3]) * 100,
    color="lightcoral",
    ls="solid",
    label="LNM (Layer 8)",
    linewidth=1,
)
ax1.plot(
    Cox_list,
    RI_anmx[:, 4] / (RI_anmx[:, 4] + RI_denit[:, 4]) * 100,
    color="mediumorchid",
    ls="solid",
    label="TNM (Layer 11)",
    linewidth=1,
)
ax1.plot(
    Cox_list,
    RI_anmx[:, 5] / (RI_anmx[:, 5] + RI_denit[:, 5]) * 100,
    color="deeppink",
    ls="solid",
    label="ODZ Bottom (Layer 14)",
    linewidth=1,
)
ax1.fill_between(np.arange(-1.92, -0.32, 0.2), 0, 100, color="gainsboro", alpha=0.5)
ax1.axvline(x=0.18, color="black", ls="--", linewidth=1)
ax1.axvline(x=0.99, color="black", ls="--", linewidth=1)
ax1.legend(loc=2, fontsize=6)
ax1.set(ylim=(20, 100), xlim=(-2.5, 1.5))
ax1.set_ylabel("Anammox %", fontsize=fs - 4)
ax1.set_xlabel("Carbon Oxidation State", fontsize=fs - 3)
ax1.tick_params(axis="x", labelsize=fs)
ax1.tick_params(axis="y", labelsize=fs)

# Reduction/Oxidation
ax2.plot(
    Cox_list,
    RI_nitox[:, 1] / (RI_anmx[:, 1] + RI_denit[:, 1]),
    color="mediumturquoise",
    ls="solid",
    label="Layer 4",
    linewidth=1,
)
ax2.plot(
    Cox_list,
    RI_nitox[:, 2] / (RI_anmx[:, 2] + RI_denit[:, 2]),
    color="springgreen",
    ls="solid",
    label="Layer 6",
    linewidth=1,
)
ax2.plot(
    Cox_list,
    RI_nitox[:, 3] / (RI_anmx[:, 3] + RI_denit[:, 3]),
    color="lightcoral",
    ls="solid",
    label="Layer 8",
    linewidth=1,
)
ax2.plot(
    Cox_list,
    RI_nitox[:, 4] / (RI_anmx[:, 4] + RI_denit[:, 4]),
    color="mediumorchid",
    ls="solid",
    label="Layer 11",
    linewidth=1,
)
ax2.plot(
    Cox_list,
    RI_nitox[:, 5] / (RI_anmx[:, 5] + RI_denit[:, 5]),
    color="deeppink",
    ls="solid",
    label="Layer 14",
    linewidth=1,
)
ax2.fill_between(np.arange(-1.92, -0.32, 0.2), 0, 16, color="gainsboro", alpha=0.5)
ax2.axvline(x=0.18, color="black", ls="--", linewidth=1)
ax2.axvline(x=0.99, color="black", ls="--", linewidth=1)
ax2.set(ylim=(0, 16), xlim=(-2.5, 1.5))
ax2.set_ylabel("Nitrite Oxidation to Reduction Ratio", fontsize=fs - 4)
ax2.set_xlabel("Carbon Oxidation State", fontsize=fs - 3)
ax2.tick_params(axis="x", labelsize=fs)
ax2.tick_params(axis="y", labelsize=fs)

# DNRN/Nitrite Oxidation
ax3.plot(
    Cox_list,
    RI_nitox[:, 1] / RI_DNRN[:, 1],
    color="mediumturquoise",
    label="Layer 4",
    linewidth=1,
)
ax3.plot(
    Cox_list,
    RI_nitox[:, 2] / RI_DNRN[:, 2],
    color="springgreen",
    label="Layer 6",
    linewidth=1,
)
ax3.plot(
    Cox_list,
    RI_nitox[:, 3] / RI_DNRN[:, 3],
    color="lightcoral",
    label="Layer 8",
    linewidth=1,
)
ax3.plot(
    Cox_list,
    RI_nitox[:, 4] / RI_DNRN[:, 4],
    color="mediumorchid",
    label="Layer 11",
    linewidth=1,
)
ax3.plot(
    Cox_list,
    RI_nitox[:, 5] / RI_DNRN[:, 5],
    color="deeppink",
    label="Layer 14",
    linewidth=1,
)
ax3.set(ylim=(0, 1.2), xlim=(-2.5, 1.5))
ax3.fill_between(np.arange(-1.92, -0.32, 0.2), 0, 1.2, color="gainsboro", alpha=0.5)
ax3.axvline(x=0.18, color="black", ls="--", linewidth=1)
ax3.axvline(x=0.99, color="black", ls="--", linewidth=1)
ax3.set_ylabel("Nitrite Oxidation to DNRN Ratio", fontsize=fs - 4)
ax3.set_xlabel("Carbon Oxidation State", fontsize=fs - 3)
ax3.tick_params(axis="x", labelsize=fs)
ax3.tick_params(axis="y", labelsize=fs)

fig.tight_layout()
plt.savefig(
    "figures/OM_variations/comparetive/ratios_layers_comp.{}".format(fig_format),
    dpi=dpi,
)
print("done")
