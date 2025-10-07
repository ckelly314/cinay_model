# Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gsw
import PyCO2SYS as pyco2

# Set Directory

# Load Datasets
CLIVAR_df = pd.read_csv(
    "../../Data/P18_cruises/CLIVAR_2007.csv", skiprows=np.arange(0, 52)
).iloc[
    1:-1
]  # Removing  headers, metadata etc.
GOSHIP_df = pd.read_csv(
    "../../Data/P18_cruises/GOSHIP_2016.csv", skiprows=np.arange(0, 76)
).iloc[1:-1]
falkor_df = pd.read_csv("falkor_clean.csv")


# Select Unflagged Data at 110 W, 14N: P, T, S, Nitrate etc.
CLIVAR = (
    CLIVAR_df.loc[
        ((CLIVAR_df["NITRIT_FLAG_W"] == 2) | (CLIVAR_df["NITRIT_FLAG_W"] == 6))
        & ((CLIVAR_df["OXYGEN_FLAG_W"] == 2) | (CLIVAR_df["OXYGEN_FLAG_W"] == 6))
        & ((CLIVAR_df["ALKALI_FLAG_W"] == 2) | (CLIVAR_df["ALKALI_FLAG_W"] == 6))
        & ((CLIVAR_df["PH_SWS_FLAG_W"] == 2) | (CLIVAR_df["PH_SWS_FLAG_W"] == 6))
        & ((CLIVAR_df["STNNBR"] == 19))
    ][
        [
            "LONGITUDE",
            "LATITUDE",
            "CTDPRS",
            "CTDTMP",
            "CTDSAL",
            "NITRIT",
            "OXYGEN",
            "PHSPHT",
            "ALKALI",
            "PH_SWS",
            "PH_TMP",
        ]
    ]
    .reset_index(drop=True)
    .astype("float")
)
GOSHIP = (
    GOSHIP_df.loc[
        ((GOSHIP_df["NITRIT_FLAG_W"] == 2) | (GOSHIP_df["NITRIT_FLAG_W"] == 6))
        & ((GOSHIP_df["OXYGEN_FLAG_W"] == 2) | (GOSHIP_df["OXYGEN_FLAG_W"] == 6))
        & ((GOSHIP_df["ALKALI_FLAG_W"] == 2) | (GOSHIP_df["ALKALI_FLAG_W"] == 6))
        & ((GOSHIP_df["PH_SWS_FLAG_W"] == 2) | (GOSHIP_df["PH_SWS_FLAG_W"] == 6))
        & ((GOSHIP_df["STNNBR"] == 21))
    ][
        [
            "LONGITUDE",
            "LATITUDE",
            "CTDPRS",
            "CTDTMP",
            "CTDSAL",
            "NITRIT",
            "OXYGEN",
            "PHSPHT",
            "ALKALI",
            "PH_SWS",
            "PH_TMP",
        ]
    ]
    .reset_index(drop=True)
    .astype("float")
)
falkor = falkor_df.loc[falkor_df["Station"] == 9][
    ["lon", "lat", "P", "sigma0", "NO2", "O2", "pH"]
].reset_index(drop=True)

# Calculate Sigma0
CLIVAR["sigma0"] = gsw.density.sigma0(
    gsw.SA_from_SP(
        CLIVAR["CTDSAL"], CLIVAR["CTDPRS"], CLIVAR["LONGITUDE"], CLIVAR["LATITUDE"]
    ),
    gsw.CT_from_t(
        gsw.SA_from_SP(
            CLIVAR["CTDSAL"], CLIVAR["CTDPRS"], CLIVAR["LONGITUDE"], CLIVAR["LATITUDE"]
        ),
        CLIVAR["CTDTMP"],
        CLIVAR["CTDPRS"],
    ),
)
GOSHIP["sigma0"] = gsw.density.sigma0(
    gsw.SA_from_SP(
        GOSHIP["CTDSAL"], GOSHIP["CTDPRS"], GOSHIP["LONGITUDE"], GOSHIP["LATITUDE"]
    ),
    gsw.CT_from_t(
        gsw.SA_from_SP(
            GOSHIP["CTDSAL"], GOSHIP["CTDPRS"], GOSHIP["LONGITUDE"], GOSHIP["LATITUDE"]
        ),
        GOSHIP["CTDTMP"],
        GOSHIP["CTDPRS"],
    ),
)


# Calculate in Situ pH
CLIVAR["pH_insitu"] = pyco2.sys(
    par1=CLIVAR["ALKALI"],
    par2=CLIVAR["PH_SWS"],
    par1_type=1,
    par2_type=3,
    salinity=CLIVAR["CTDSAL"],
    temperature=CLIVAR["PH_TMP"],
    pressure=0,
    temperature_out=CLIVAR["CTDTMP"],
    pressure_out=CLIVAR["CTDPRS"],
    total_phosphate=CLIVAR["PHSPHT"],
    opt_pH_scale=2,
    opt_k_carbonic=15,
)["pH_out"]
GOSHIP["pH_insitu"] = pyco2.sys(
    par1=GOSHIP["ALKALI"],
    par2=GOSHIP["PH_SWS"],
    par1_type=1,
    par2_type=3,
    salinity=GOSHIP["CTDSAL"],
    temperature=GOSHIP["PH_TMP"],
    pressure=0,
    temperature_out=GOSHIP["CTDTMP"],
    pressure_out=GOSHIP["CTDPRS"],
    total_phosphate=GOSHIP["PHSPHT"],
    opt_pH_scale=2,
    opt_k_carbonic=15,
)["pH_out"]


# Falkor Station 9
falkor.sort_values(by="sigma0", inplace=True, ignore_index=True)
st9 = falkor.loc[(falkor["sigma0"] >= 25) & (falkor["sigma0"] <= 27.4)]

# Plot Vertical NO2- Profiles (vs. Sigma0 and depth)
fig, (ax2, ax1, ax3) = plt.subplots(1, 3, figsize=(10, 6))

ax1.plot(
    st9["NO2"],
    st9["sigma0"],
    color="lightcoral",
    marker="o",
    ms=3,
    linewidth=1,
    label="FK180624 St 9, 2018",
)
ax1.plot(
    CLIVAR["NITRIT"],
    CLIVAR["sigma0"],
    color="mediumturquoise",
    marker="v",
    ms=3,
    linewidth=1,
    label="CLIVAR St 19, 2007",
)
ax1.plot(
    GOSHIP["NITRIT"],
    GOSHIP["sigma0"],
    color="mediumorchid",
    marker="*",
    ms=3,
    linewidth=1,
    label="GOSHIP St 21, 2016",
)
ax1.set(ylim=(25, 27.4), xlim=(0, 3))
ax1.xaxis.set_label_position("top")
ax1.xaxis.set_ticks_position("top")
ax1.invert_yaxis()
# ax1.set_ylabel('$\sigma$$_{\\theta}$ ($kgm^{-3}$)', fontsize = 13)
ax1.set_xlabel("[$NO_2^-$] ($\mu$$mol$ $kg^{-1})$", fontsize=15)
ax1.axes.yaxis.set_ticklabels([])
ax1.set_xticks(ticks=np.array((0.0, 1.0, 2.0, 3.0)))
ax1.tick_params(axis="x", labelsize=15)
ax1.legend()

ax2.plot(
    st9["O2"] - 0.5,
    st9["sigma0"],
    color="lightcoral",
    marker="o",
    ms=3,
    linewidth=1,
    label="Falkor St 9, 2018",
)
ax2.plot(
    CLIVAR["OXYGEN"] - 3.2,
    CLIVAR["sigma0"],
    color="mediumturquoise",
    marker="v",
    ms=3,
    linewidth=1,
    label="CLIVAR St 19, 2007",
)
ax2.plot(
    GOSHIP["OXYGEN"] - 2.5,
    GOSHIP["sigma0"],
    color="mediumorchid",
    marker="*",
    ms=3,
    linewidth=1,
    label="GOSHIP St 21, 2016",
)
ax2.set(ylim=(25, 27.4), xlim=(-0.5, 15))
ax2.xaxis.set_label_position("top")
ax2.xaxis.set_ticks_position("top")
ax2.invert_yaxis()
ax2.set_ylabel("$\sigma$$_{\\theta}$ ($kg$ $m^{-3}$)", fontsize=15)
ax2.set_xlabel("$\\tilde{O}$$_2$ ($\mu$$mol$ $kg^{-1})$", fontsize=15)
ax2.tick_params(axis="x", labelsize=15)
ax2.tick_params(axis="y", labelsize=15)
ax2.set_xticks(ticks=np.array((0.0, 5.0, 10.0, 15.0)))


ax3.plot(
    st9["pH"],
    st9["sigma0"],
    color="lightcoral",
    marker="o",
    ms=3,
    linewidth=1,
    label="Falkor St 9, 2018",
)
ax3.plot(
    CLIVAR["pH_insitu"],
    CLIVAR["sigma0"],
    color="mediumturquoise",
    marker="v",
    ms=3,
    linewidth=1,
    label="CLIVAR St 19, 2007",
)
ax3.plot(
    GOSHIP["pH_insitu"],
    GOSHIP["sigma0"],
    color="mediumorchid",
    marker="*",
    ms=3,
    linewidth=1,
    label="GOSHIP St 21, 2016",
)
ax3.set(ylim=(25, 27.4), xlim=(7.52, 7.64))
ax3.xaxis.set_label_position("top")
ax3.xaxis.set_ticks_position("top")
ax3.invert_yaxis()
# ax3.set_ylabel('$\sigma$$_{\\theta}$ ($kgm^{-3}$)')
ax3.axes.yaxis.set_ticklabels([])
ax3.set_xlabel("pH", fontsize=15)
ax3.tick_params(axis="x", labelsize=15)
ax3.set_xticks(ticks=np.arange(7.52, 7.65, 0.04))

fig.tight_layout()
plt.savefig("figures/station_plots/multiple_P18_profiles.png", dpi=500)
