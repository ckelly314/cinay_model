"""
File: pH_and_reactions.py
----------------------------

Calculates the pH increase per organic carbon remineralized by denitrification
for various organic matter compositions with varying carbon oxidation states
for a range of DIC and TA values observed in the ETNP ODZ.
Plots delta_pH versus carbon oxidation state.
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PyCO2SYS as pyco2
from matplotlib.patches import Rectangle


# Load Data
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="Station")

# Define Variable
TA_range = np.linspace(2200, 2450, 200)  # umol/kg
DIC_range = np.linspace(2000, 2400, 200)  # umol/kg
pH_range = np.linspace(7, 8, 1000)
DIC_box, TA_box = np.meshgrid(DIC_range, TA_range)
T = 11  # degrees C
S = 35  # psu
P = 300  # dbar

# Reaction Distribution
ri_denit = 1
dTA_box = 1.096 * ri_denit
dDIC_box = ri_denit
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

# Contour Map
fig, ax = plt.subplots(figsize=(9, 7))
plt.contourf(DIC_box, TA_box, pH_new - pH_baseline, vmin=0, vmax=0.0004)
ax.add_patch(Rectangle((2250, 2290), 50, 25, fc="none", color="red", linewidth=2))
cb = plt.colorbar()
plt.xlabel("DIC ($\mu$$mol$ $kg^{-1}$)", fontsize=16)
plt.ylabel("Total Alkalinity ($\mu$$mol$ $kg^{-1}$)", fontsize=16)
cb.set_label("$\Delta$pH per $C_{org}$ remineralized", fontsize=16)
plt.tick_params(axis="x", labelsize=16)
plt.tick_params(axis="y", labelsize=16)
cb.ax.tick_params(labelsize=16)
plt.savefig("figures/delta_pH_contour_denit_new.pdf")
plt.show()
print("Done")
