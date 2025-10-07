"""
File: pH_and_reactions_Cox.py
-----------------------------

Calculates the pH increase per organic carbon remineralized by denitrification
for the estimated organic matter composition based on FK180624 data.
Maps delta_pH for a range of DIC and TA values.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PyCO2SYS as pyco2


# Load Data
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="Station")
C_ox = np.arange(-4, 4.5, 0.5)
data_path = "output/OM_variations/{}/{}"
fpath = list()
for i in C_ox:
    if i < 0:
        fpath.append("cox_neg{}".format(round(abs(i), 2)))
    elif i == 0:
        fpath.append("cox_zero")
    elif i > 0:
        fpath.append("cox_pos{}".format(round(i, 1)))
RI_anmx = np.zeros((len(C_ox), 5))
RI_denit = np.zeros((len(C_ox), 5))
for i in np.arange(0, len(C_ox)):
    nit_sys = np.array(
        pd.read_csv(data_path.format(fpath[i], "relative_importances_nitrite_mean.csv"))
    )
    RI_anmx[i, :] = np.array(
        [
            np.mean(nit_sys[:, 0]),
            nit_sys[3, 0],
            nit_sys[6, 0],
            nit_sys[10, 0],
            nit_sys[13, 0],
        ]
    )
    RI_denit[i, :] = np.array(
        [
            np.mean(nit_sys[:, 1]),
            nit_sys[3, 1],
            nit_sys[6, 1],
            nit_sys[10, 1],
            nit_sys[13, 1],
        ]
    )

# Define Variable
TA_range = np.linspace(2290, 2315, 100)  # umol/kg
DIC_range = np.linspace(2250, 2300, 100)  # umol/kg
DIC_box, TA_box = np.meshgrid(DIC_range, TA_range)
T = 11  # degrees C
S = 35  # psu
P = 300  # dbar

# Carbon Oxidation State Variations (only denitrification)
dpH_denit = np.zeros((100 * 100, np.size(C_ox)))
for i in np.arange(0, np.size(C_ox)):
    C_ox_iter = C_ox[i]
    charge = 0
    # Stoichiometric Ratios
    rN = 11.5
    rC = 112.5
    twoO_minus_H = C_ox_iter * rC - charge - 3 * rN + 5 * 1
    dTA_box = ((4 / 3 * rC - 1 / 3 * twoO_minus_H - rN + 5 / 3) + rN - 1) / rC
    dDIC_box = 1

    # Calculate pH using CO2SYS
    pH_baseline = pyco2.sys(
        par1=TA_box,
        par2=DIC_box,
        par1_type=1,
        par2_type=2,
        salinity=S,
        temperature=T,
        pressure=P,
        total_phosphate=2.5,
    )["pH"]
    pH_new = pyco2.sys(
        par1=TA_box + dTA_box,
        par2=DIC_box + dDIC_box,
        par1_type=1,
        par2_type=2,
        salinity=S,
        temperature=T,
        pressure=P,
        total_phosphate=2.5,
    )["pH"]
    dpH_iter = pH_new - pH_baseline
    dpH_denit[:, i] = dpH_iter.flatten()
df_denit = pd.DataFrame(dpH_denit, columns=C_ox)
mean_denit = df_denit.mean()
sd_denit = df_denit.std()

# Plot
fig, ax = plt.subplots(figsize=(9, 7))
plt.errorbar(
    C_ox,
    mean_denit,
    yerr=sd_denit,
    linewidth=1.5,
    marker="o",
    markerfacecolor="white",
    ms=7,
    capsize=5,
    ecolor="black",
    barsabove=True,
    color="mediumturquoise",
)
plt.axhline(y=0.00, color="black", ls="-", linewidth=1)
plt.fill_between(
    np.arange(-1.92, -0.32, 0.2), -0.04, 0.06, color="gainsboro", alpha=0.5
)
plt.axvline(x=0.18, color="black", ls="--", linewidth=1)
plt.axvline(x=0.99, color="black", ls="--", linewidth=1)

plt.xlabel("Carbon Oxidation State", fontsize=16)
plt.ylabel("$\Delta$pH per $C_{org}$ remineralized", fontsize=16)
plt.tick_params(axis="x", labelsize=16)
plt.tick_params(axis="y", labelsize=16)
plt.xlim((-2.5, 1.5))
plt.ylim((-0.004, 0.006))
fig.tight_layout()
plt.savefig("figures/delta_pH_Cox_falkor.pdf")  # , dpi = 500)
