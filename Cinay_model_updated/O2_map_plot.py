# Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cf
import xarray as xr


# Load O2 Map Data
O2 = np.loadtxt("output/O2_map/O2_nonan0.csv", delimiter=",")
lon = np.loadtxt("output/O2_map/O2_lon.csv", delimiter=",")
lat = np.loadtxt("output/O2_map/O2_lat.csv", delimiter=",")

# Create Xarray File
O2_data = xr.DataArray(O2, coords=[lon, lat], dims=["lon", "lat"]).T
O2_interp = O2_data.interp(
    lon=np.arange(np.min(lon), np.max(lon) + 0.1, 0.1),
    lat=np.arange(np.min(lat), np.max(lat) + 0.1, 0.1),
)

# Large Map O2
fig = plt.figure(figsize=(12, 4))
plt.rc("font", size=20)
plt.rcParams["font.family"] = "arial"
ax = plt.subplot(projection=ccrs.PlateCarree())
m = O2_interp.plot.contourf(
    levels=np.arange(0, 54, 10),
    cmap="cividis",
    yticks=np.arange(0, 50, 10),
    xticks=np.arange(-140, -60, 20),
    xlim=(-140, -80),
    ylim=(-5, 35),
    ax=ax,
    add_colorbar=False,
)
ax.set(xlabel="Longitude (˚E)", ylabel="Latitude (˚N)")
cb = plt.colorbar(m, ax=ax)
cb.set_label("[$O_2$] ($\mu$$mol$ $kg^{-1}$)")
ax.add_feature(cf.LAND, zorder=20, edgecolor="black", color="white")

ax.add_patch(
    Rectangle((-120, 10), 20, 10, fc="none", color="white", linewidth=2.7, ls="--")
)
ax.add_patch(Rectangle((-98, 18), 17.5, 17, fc="gray", color="gray"))
ax.add_patch(Rectangle((-90, 16), 9.5, 2, fc="gray", color="gray"))
ax.add_patch(Rectangle((-84, 9.2), 3.5, 6.7, fc="gray", color="gray"))

fig.tight_layout()
plt.savefig("figures/track/o2_map_large.png", dpi=500)

# Small O2 Map with Stations
fig = plt.figure(figsize=(12, 4))
plt.rcParams["font.family"] = "arial"
plt.rc("font", size=20)
ax = plt.subplot(projection=ccrs.PlateCarree())
m = O2_interp.plot.contourf(
    levels=np.arange(0, 54, 10),
    cmap="cividis",
    yticks=np.arange(12, 20, 2),
    xticks=np.arange(-120, -6, 4),
    xlim=(-120, -100),
    ylim=(11, 19),
    ax=ax,
    add_colorbar=False,
)
ax.set(xlabel="Longitude (˚E)", ylabel="Latitude (˚N)")
cb = plt.colorbar(m, ax=ax)
cb.set_label("[$O_2$] ($\mu$$mol$ $kg^{-1}$)")
ax.add_feature(cf.LAND, zorder=20, edgecolor="black", color="white")

# Add Stations
falkor = pd.read_csv("falkor_clean.csv").sort_values(by="P")
station = np.array(falkor["Station"])
lon_st = np.array(falkor["lon"])
lat_st = np.array(falkor["lat"])
station_vec = np.arange(1, 20, 1)
lon_vec = np.ones(np.shape(station_vec))
lat_vec = np.ones(np.shape(station_vec))
for i in np.arange(1, 20, 1):
    idx = np.where(station == i)
    lon_vec[i - 1] = np.mean(lon_st[idx])
    lat_vec[i - 1] = np.mean(lat_st[idx])
ax.scatter(lon_vec, lat_vec, c="orange", s=50)
for i in np.arange(1, 20):
    ax.text(
        lon_vec[i - 1],
        lat_vec[i - 1] - 0.7,
        "{}".format(i),
        color="orange",
        fontsize=18,
        horizontalalignment="center",
    )
fig.tight_layout()
plt.savefig("figures/track/o2_map_track.png", dpi=500)
print("done")
