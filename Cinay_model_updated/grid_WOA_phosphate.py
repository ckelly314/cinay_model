import xarray as xr
import pandas as pd
import numpy as np

def round_off(series):
	return round(series * 2) / 2

def customround(series, base=5):
	return base * round(series/base)

f5906484 = pd.read_csv("5906484qcno2.txt",skiprows = 70, sep='\t',
    parse_dates = ['mon/day/yr'])
# add 'month' column so we can match gridded WOA data
f5906484['month'] = f5906484['mon/day/yr'].dt.strftime('%m')
f5906484['month'] = f5906484['month'].astype(int)

data = f5906484

data['depth'] = customround(data['Depth[m]'])
data['lon'] = round(data['Lon [°E]']) + 0.5 #round_off(data['Lon [°E]'])
data['lon'] = data['lon'] - 360.0
data['lat'] = round(data['Lat [°N]']) + 0.5 #round_off(data['Lat [°N]'])


def process_WOA(month):
	if i<10:
		link = f'https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/phosphate/all/1.00/woa18_all_p0{month}_01.nc'
	elif i>=10:
		link = f'https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/phosphate/all/1.00/woa18_all_p{month}_01.nc'
	data = xr.open_dataset(link, decode_times = False)
	data = data.sel(lat = slice(16,19),lon = slice(-110, -106))
	data = data.to_dataframe(dim_order = ['lat', 'nbounds', 'lon', 'depth', 'time'])
	data = data.reset_index()
	data = data[['lat','lon','depth','p_an']].groupby(['lat','lon','depth']).mean()
	data['month'] = month
	return data

jan = xr.open_dataset("https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/phosphate/all/1.00/woa18_all_p01_01.nc", decode_times = False)
jan = jan.sel(lat = slice(16,19),lon = slice(-110, -106))
jan = jan.to_dataframe(dim_order = ['lat', 'nbounds', 'lon', 'depth', 'time'])
jan = jan.reset_index()
jan['month'] = 1
#jan = jan[['lat','lon','depth','p_an']].groupby(['lat','lon','depth']).mean()

data = data.set_index(['month','lat','lon','depth'])

WOA = pd.DataFrame()

for i in range(1,13):
	print(i)
	monthlyWOA = process_WOA(i)
	WOA = pd.concat([WOA, monthlyWOA])

WOA = WOA.reset_index()
WOA = WOA.set_index(['month','lat','lon','depth'])

data = data.join(WOA)

print(data['p_an'].head())
