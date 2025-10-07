"""
File: calcium_carbonate_saturation.py
-------------------------------------

Fits observed pH values from FK180624 data and creates a hypothetical pH profile.
Calculates omega for using these pH profiles and in-situ temperature, alkalinity, salinity, etc. profiles.
Plots both pH and omega as a function of depth. 
"""

## Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PyCO2SYS as pyco2
import seaborn as sns
import scipy as sc
from scipy.signal import savgol_filter

# Figure Settings
dpi = 500
ff = "arial"
fs = 12
fig_format = "pdf"

# Load Data and Define Tracers
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="P")
st = np.array((2))
T = np.array(falkor["T"])  # degrees C
S = np.array(falkor["S"])  # psu
P = np.array(falkor["P"])  # dbar
sigma0 = np.array(falkor["sigma0"])  # kg/m3
TA = np.array(falkor["TA"])  # umol/kg
DIP = np.array(falkor["DIP"])  # umol/kg
NH4 = np.array(falkor["NH4"])  # umol/kg
# silicate = np.array(falkor['silicate']) # umol/kg
pH = np.array(falkor["pH"])
O2 = np.array(falkor["O2"])  # umol/kg

# Create pH Profiles
N = 100000
p0 = np.array([1.14, 0.02, -0.11, 108.83, 26.72, 7.55])


def pH_fit(depth, a, b, c, d, e, f):
    pH_est = (
        a * np.exp(-b * (depth)) + c * np.exp(-((depth - d) ** 2) / (2 * e**2))
    ) + f
    return pH_est


def pH_fit2(depth, a, b, c, f):
    pH_est = a * (np.exp(-b * (-depth - c))) + f
    return pH_est


depth = P[np.where((P >= 50) & (P <= 900))]
pH_obs = pH[np.where((P >= 50) & (P <= 900))]
popt = sc.optimize.curve_fit(pH_fit, depth, pH_obs, p0=p0)
depth_vec = np.linspace(50, 900, N)
pH_vec1 = pH_fit(depth_vec, *np.array(popt[0]))
pH_vec2 = pH_fit2(depth_vec, 6, -0.05, -1, 7.556)


# Smooth and Interpolate CO2SYS inputs
T_smooth = savgol_filter(T[np.where((P >= 50) & (P <= 900))], 101, 3)
T_vec = sc.interpolate.interp1d(
    depth, T_smooth, kind="linear", fill_value="extrapolate"
)(depth_vec)
S_smooth = savgol_filter(S[np.where((P >= 50) & (P <= 900))], 101, 3)
S_vec = sc.interpolate.interp1d(
    depth, S_smooth, kind="linear", fill_value="extrapolate"
)(depth_vec)
DIP_smooth = savgol_filter(DIP[np.where((P >= 50) & (P <= 900))], 101, 3)
DIP_vec = sc.interpolate.interp1d(
    depth, DIP_smooth, kind="linear", fill_value="extrapolate"
)(depth_vec)
NH4_smooth = savgol_filter(NH4[np.where((P >= 50) & (P <= 900))], 101, 3)
NH4_vec = sc.interpolate.interp1d(
    depth, NH4_smooth, kind="linear", fill_value="extrapolate"
)(depth_vec)
TA_smooth = savgol_filter(TA[np.where((P >= 50) & (P <= 900))], 101, 3)
TA_vec = sc.interpolate.interp1d(
    depth, TA_smooth, kind="linear", fill_value="extrapolate"
)(depth_vec)
# Si_smooth = savgol_filter(silicate[np.where((P >= 50) & (P <= 900))], 101,3)
# Si_vec = sc.interpolate.interp1d(depth, Si_smooth, kind ='linear', fill_value = 'extrapolate')(depth_vec)

# Run Calculate Omegas with CO2SYS
Z1 = pyco2.sys(
    par1=TA_vec,
    par2=pH_vec1,
    par1_type=1,
    par2_type=3,
    salinity=S_vec,
    temperature=T_vec,
    # pressure = depth_vec, total_silicate= Si_vec, total_phosphate = DIP_vec, total_ammonia = NH4_vec, opt_pH_scale = 1)
    pressure=depth_vec,
    total_phosphate=DIP_vec,
    total_ammonia=NH4_vec,
    opt_pH_scale=1,
)
omegaA_vec1 = np.array(Z1["saturation_aragonite"], dtype="float")
omegaC_vec1 = np.array(Z1["saturation_calcite"], dtype="float")

Z2 = pyco2.sys(
    par1=TA_vec,
    par2=pH_vec2,
    par1_type=1,
    par2_type=3,
    salinity=S_vec,
    temperature=T_vec,
    # pressure = depth_vec, total_silicate= Si_vec, total_phosphate = DIP_vec, total_ammonia = NH4_vec, opt_pH_scale = 1)
    pressure=depth_vec,
    total_phosphate=DIP_vec,
    total_ammonia=NH4_vec,
    opt_pH_scale=1,
)
omegaA_vec2 = np.array(Z2["saturation_aragonite"], dtype="float")
omegaC_vec2 = np.array(Z2["saturation_calcite"], dtype="float")

# Plot Different Profiles
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 7))
sns.scatterplot(x=pH, y=P, s=7, ax=ax1)
ax1.plot(pH_vec1, depth_vec, color="black", ls="-", linewidth=1.5)
ax1.plot(pH_vec2, depth_vec, color="black", ls="--", linewidth=1.5)
ax1.set_ylim((0, 600))
ax1.set_xlim((7.54, 7.65))
ax1.invert_yaxis()
ax1.set_xlabel("pH", fontsize=fs)
ax1.set_ylabel("Depth (dbar)", fontsize=fs)
ax1.xaxis.set_label_position("top")
ax1.xaxis.set_ticks_position("top")
ax1.tick_params(axis="x", labelsize=fs)
ax1.tick_params(axis="y", labelsize=fs)

ax2.plot(
    omegaA_vec1,
    depth_vec,
    color="mediumturquoise",
    ls="-",
    linewidth=1.5,
    label="Aragonite",
)
ax2.plot(omegaA_vec2, depth_vec, color="mediumturquoise", ls="--", linewidth=1.5)
ax2.plot(
    omegaC_vec1, depth_vec, color="lightcoral", ls="-", linewidth=1.5, label="Calcite"
)
ax2.plot(omegaC_vec2, depth_vec, color="lightcoral", ls="--", linewidth=1.5)
ax2.axvline(x=1, c="black", ls=":")
ax2.set_ylim((0, 600))
ax2.set_xlim((0.5, 1.6))
ax2.invert_yaxis()
ax2.legend()
ax2.set_xlabel("$\Omega$", fontsize=fs)
# ax2.set_ylabel('Depth (dbar)', fontsize = fs)
ax2.xaxis.set_label_position("top")
ax2.xaxis.set_ticks_position("top")
ax2.tick_params(axis="x", labelsize=fs)
ax2.axes.yaxis.set_ticklabels([])


fig1.tight_layout(pad=2.5)
plt.savefig("figures/model_pH_profile.{}".format(fig_format))  # , dpi = dpi)
