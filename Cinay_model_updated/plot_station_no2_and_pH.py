## Import Libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Set Directory
img_format = ".png"
fs = 10
ff = "arial"

# Load Data and Define Variables
falkor = pd.read_csv("falkor_clean.csv")
NO2 = np.array(falkor["NO2"])
st = np.array(falkor["Station"])
sigma0 = np.array(falkor["sigma0"])
pH = np.array(falkor["pH"])
O2 = np.array(falkor["O2"])
depth = np.array(falkor["P"])

# Calculate the Maximum Values for Each Station
max_pH = np.ones(np.shape(st))
max_no2 = np.ones(np.shape(st))
max_no2_unique = np.zeros(np.shape(np.arange(1, 20)))
max_pH_unique = np.zeros(np.shape(np.arange(1, 20)))
for i in np.arange(1, 20):
    idx = np.where((st == i) & (sigma0 > 25.4) & (sigma0 < 27.00))
    max_no2[idx] = np.round(np.max(NO2[idx]), 2)
    max_no2_unique[i - 1] = np.round(np.max(NO2[idx]), 2)
    max_pH[idx] = np.round(np.max(pH[idx]), 4)
    max_pH_unique[i - 1] = np.round(np.max(pH[idx]), 4)
df = pd.DataFrame(
    {
        "station": st,
        "sigma0": sigma0,
        "NO2": NO2,
        "pH": pH,
        "max_no2": max_no2,
        "max_pH": max_pH,
        "O2": O2,
        "depth": depth,
    }
).loc[np.where((sigma0 > 25.4) & (sigma0 < 27))]
df_full = pd.DataFrame(
    {"station": st, "sigma0": sigma0, "NO2": NO2, "pH": pH, "depth": depth}
).loc[np.where((sigma0 > 21) & (sigma0 < 28))]
df.sort_values(by=["station", "sigma0"], inplace=True)
print(df_full.describe())
# Create the color vector
palette = cm.get_cmap("Blues", 25)
colors = palette(np.linspace(0, 1, 25))[2:22]
df_color = pd.DataFrame({"station": np.arange(1, 20), "max_no2": max_no2_unique})
df_color.sort_values(by=["max_no2"], ascending=True, inplace=True, ignore_index=True)
index = df_color.index[df_color["station"] == 17].tolist()

# Plot Nitrite Profile
fig1, ax = plt.subplots(1, 19, figsize=(10, 2.5))
for i in np.arange(0, 19):
    data = df.iloc[np.where(np.array(df["station"]) == int(i + 1))]
    fill_color = colors[df_color.index[df_color["station"] == i + 1].tolist()]
    ax[-(i + 1)].plot(
        data["NO2"], data["sigma0"], linewidth=1.0, color="black"
    )  # Change outline properties
    ax[-(i + 1)].fill_betweenx(data["sigma0"], data["NO2"], 0, color=fill_color)
    ax[-(i + 1)].axvline(x=0, linewidth=1.5, color="black")  # Change Outline Properties
    ax[-(i + 1)].set(xlim=(0, 2.5), ylim=(25.4, 27))
    ax[-(i + 1)].invert_yaxis()
    ax[-(i + 1)].set_zorder(i)
    ax[1].set_zorder(18)
    ax[0].set_zorder(17)
    if i == 18:
        ax[-(i + 1)].xaxis.set_label_position("top")
        ax[-(i + 1)].xaxis.set_ticks_position("top")
        ax[-(i + 1)].spines["right"].set_visible(False)
        ax[-(i + 1)].spines["bottom"].set_visible(False)
        ax[-(i + 1)].set_xlabel("[$NO_2^-$] ($\mu$$mol$ $kg^{-1}$)", fontsize=fs)
        ax[-(i + 1)].set_ylabel("$\sigma$$_{\\theta}$ ($kg$ $m^{-3}$)", fontsize=fs)
        ax[-(i + 1)].set_yticks(ticks=np.arange(25.4, 27.2, 0.4))
        ax[-(i + 1)].set_xticks(ticks=np.array((0, 1.0, 2.5)))
        ax[-(i + 1)].set(xlim=(0, 2.5), ylim=(25.4, 27))
        ax[-(i + 1)].invert_yaxis()
        ax[-(i + 1)].tick_params(axis="x", labelsize=fs)
        ax[-(i + 1)].tick_params(axis="y", labelsize=fs)
        ax[-(i + 1)].set_title(
            " St. {}".format(i + 1), y=-0.01, fontsize=fs, loc="left"
        )
    else:
        ax[-(i + 1)].tick_params(
            left=False, bottom=False, labelleft=False, labelbottom=False
        )
        ax[-(i + 1)].axis("off")
        ax[-(i + 1)].set_title("    {}".format(i + 1), y=-0.01, fontsize=fs, loc="left")


fig1.tight_layout(pad=1.08, w_pad=-3)
plt.savefig(
    "figures/station_plots/nitrite_stations_overlapped{}".format(img_format), dpi=500
)

# Plot pH  Profile
palette = cm.get_cmap("Greens", 25)
colors = palette(np.linspace(0, 1, 25))[2:22]
df_color = pd.DataFrame({"station": np.arange(1, 20), "max_pH": max_no2_unique})
df_color.sort_values(by=["max_pH"], ascending=True, inplace=True, ignore_index=True)
index = df_color.index[df_color["station"] == 17].tolist()

fig2, ax = plt.subplots(1, 19, figsize=(10, 2.5))

for i in np.arange(0, 19):
    data = df.iloc[np.where(np.array(df["station"]) == int(i + 1))]
    fill_color = colors[df_color.index[df_color["station"] == i + 1].tolist()]
    ax[-(i + 1)].plot(data["pH"], data["sigma0"], linewidth=1, color="black")
    ax[-(i + 1)].fill_betweenx(data["sigma0"], data["pH"], 0, color=fill_color)
    ax[-(i + 1)].axvline(x=7.548, linewidth=1.5, color="black")
    ax[-(i + 1)].set(xlim=(7.548, 7.62), ylim=(25.4, 27))
    ax[-(i + 1)].invert_yaxis()
    ax[-(i + 1)].set_zorder(i)
    ax[1].set_zorder(18)
    ax[0].set_zorder(17)

    if i == 18:
        ax[-(i + 1)].xaxis.set_label_position("top")
        ax[-(i + 1)].xaxis.set_ticks_position("top")
        ax[-(i + 1)].spines["right"].set_visible(False)
        ax[-(i + 1)].spines["bottom"].set_visible(False)
        ax[-(i + 1)].set_xlabel("pH", fontsize=fs)
        ax[-(i + 1)].set_ylabel("$\sigma$$_{\\theta}$ ($kg$ $m^{-3}$)", fontsize=fs)
        ax[-(i + 1)].set_yticks(ticks=np.arange(25.4, 27.2, 0.4))
        ax[-(i + 1)].set_xticks(ticks=np.array((7.55, 7.60)))
        ax[-(i + 1)].set(xlim=(7.548, 7.62), ylim=(25.4, 27))
        ax[-(i + 1)].invert_yaxis()
        ax[-(i + 1)].tick_params(axis="x", labelsize=fs)
        ax[-(i + 1)].tick_params(axis="y", labelsize=fs)
        ax[-(i + 1)].set_title(
            " St. {}".format(i + 1), y=-0.01, fontsize=fs, loc="left"
        )

    else:
        ax[-(i + 1)].tick_params(
            left=False, bottom=False, labelleft=False, labelbottom=False
        )
        ax[-(i + 1)].axis("off")
        ax[-(i + 1)].set_title("    {}".format(i + 1), y=-0.01, fontsize=fs, loc="left")


fig2.tight_layout(pad=1.08, w_pad=-3)
plt.savefig("figures/station_plots/pH_sigma0_overlapped{}".format(img_format), dpi=500)


## Depth Versions ##
df.sort_values(by=["station", "depth"], inplace=True)
# Plot Nitrite Profile
palette = cm.get_cmap("Blues", 25)
colors = palette(np.linspace(0, 1, 25))[2:22]
df_color = pd.DataFrame({"station": np.arange(1, 20), "max_no2": max_no2_unique})
df_color.sort_values(by=["max_no2"], ascending=True, inplace=True, ignore_index=True)
index = df_color.index[df_color["station"] == 17].tolist()

fig3, ax = plt.subplots(1, 19, figsize=(10, 2.5))
for i in np.arange(0, 19):
    data = df.iloc[np.where(np.array(df["station"]) == int(i + 1))]
    fill_color = colors[df_color.index[df_color["station"] == i + 1].tolist()]
    ax[-(i + 1)].plot(data["NO2"], data["depth"], linewidth=1.0, color="black")
    ax[-(i + 1)].fill_betweenx(data["depth"], data["NO2"], 0, color=fill_color)
    ax[-(i + 1)].axvline(x=0, linewidth=1.5, color="black")
    ax[-(i + 1)].set(xlim=(0, 2.5), ylim=(90, 450))
    ax[-(i + 1)].invert_yaxis()
    ax[-(i + 1)].set_zorder(i)
    ax[1].set_zorder(18)
    ax[0].set_zorder(17)
    if i == 18:
        ax[-(i + 1)].xaxis.set_label_position("top")
        ax[-(i + 1)].xaxis.set_ticks_position("top")
        ax[-(i + 1)].spines["right"].set_visible(False)
        ax[-(i + 1)].spines["bottom"].set_visible(False)
        ax[-(i + 1)].set_xlabel("[$NO_2^-$] ($\mu$$mol$ $kg^{-1}$)", fontsize=fs)
        ax[-(i + 1)].set_ylabel("Depth (dbar)", fontsize=fs)
        ax[-(i + 1)].set_yticks(ticks=np.arange(90, 451, 90))
        ax[-(i + 1)].set_xticks(ticks=np.array((0, 1.0, 2.5)))
        ax[-(i + 1)].set(xlim=(0, 2.5), ylim=(90, 450))
        ax[-(i + 1)].invert_yaxis()
        ax[-(i + 1)].tick_params(axis="x", labelsize=fs)
        ax[-(i + 1)].tick_params(axis="y", labelsize=fs)
        ax[-(i + 1)].set_title(
            " St. {}".format(i + 1), y=-0.01, fontsize=fs, loc="left"
        )
    else:
        ax[-(i + 1)].tick_params(
            left=False, bottom=False, labelleft=False, labelbottom=False
        )
        ax[-(i + 1)].axis("off")
        ax[-(i + 1)].set_title("    {}".format(i + 1), y=-0.01, fontsize=fs, loc="left")


fig3.tight_layout(pad=1.08, w_pad=-3)
plt.savefig(
    "figures/station_plots/nitrite_stations_depth_overlapped{}".format(img_format),
    dpi=500,
)

# pH
palette = cm.get_cmap("Greens", 25)
colors = palette(np.linspace(0, 1, 25))[2:22]
df_color = pd.DataFrame({"station": np.arange(1, 20), "max_pH": max_no2_unique})
df_color.sort_values(by=["max_pH"], ascending=True, inplace=True, ignore_index=True)
index = df_color.index[df_color["station"] == 17].tolist()

fig4, ax = plt.subplots(1, 19, figsize=(10, 2.5))
for i in np.arange(0, 19):
    data = df.iloc[np.where(np.array(df["station"]) == int(i + 1))]
    fill_color = colors[df_color.index[df_color["station"] == i + 1].tolist()]
    ax[-(i + 1)].plot(data["pH"], data["depth"], linewidth=1.0, color="black")
    ax[-(i + 1)].fill_betweenx(data["depth"], data["pH"], 0, color=fill_color)
    ax[-(i + 1)].axvline(x=7.548, linewidth=1.5, color="black")
    ax[-(i + 1)].set(xlim=(7.548, 7.62), ylim=(90, 450))
    ax[-(i + 1)].invert_yaxis()
    ax[-(i + 1)].set_zorder(i)
    ax[1].set_zorder(18)
    ax[0].set_zorder(17)

    if i == 18:
        ax[-(i + 1)].xaxis.set_label_position("top")
        ax[-(i + 1)].xaxis.set_ticks_position("top")
        ax[-(i + 1)].spines["right"].set_visible(False)
        ax[-(i + 1)].spines["bottom"].set_visible(False)
        ax[-(i + 1)].set_xlabel("pH", fontsize=fs)
        ax[-(i + 1)].set_ylabel("Depth (dbar)", fontsize=fs)
        ax[-(i + 1)].set_yticks(ticks=np.arange(90, 451, 90))
        ax[-(i + 1)].set_xticks(ticks=np.array((7.55, 7.60)))
        ax[-(i + 1)].set(xlim=(7.548, 7.62), ylim=(90, 450))
        ax[-(i + 1)].invert_yaxis()
        ax[-(i + 1)].tick_params(axis="x", labelsize=fs)
        ax[-(i + 1)].tick_params(axis="y", labelsize=fs)
        # ax[-(i+1)].set_title(' St. {}'.format(i+1), y = -0.01, fontsize = fs, loc = 'left')
    else:
        ax[-(i + 1)].tick_params(
            left=False, bottom=False, labelleft=False, labelbottom=False
        )
        ax[-(i + 1)].axis("off")
        # ax[-(i+1)].set_title('    {}'.format(i+1), y = -0.01, fontsize = fs, loc = 'left')


fig4.tight_layout(pad=1.08, w_pad=-3)
plt.savefig("figures/station_plots/pH_depth_overlapped{}".format(img_format), dpi=500)


print("done")
