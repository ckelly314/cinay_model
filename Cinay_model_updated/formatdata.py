import pandas as pd
from datetime import datetime

data = pd.read_csv("5906484qcno2_updated.txt",skiprows = 79, sep='\t',
    parse_dates = ['mon/day/yr'], na_values = [-10000000000.0])

data.rename(columns={
    'Lon [°E]': 'lon',
    'Lat [°N]': 'lat',
    'Pressure[dbar]':'pres',
    'Temperature[°C]': 'temp',
    'Salinity[pss]': 'psal',
    'Oxygen[µmol/kg]': 'doxy'}, inplace=True)


data['gtime'] = pd.to_datetime(
    data['mon/day/yr'].astype(str) + ' ' + data['hh:mm'].astype(str),
    format="%Y-%m-%d %H:%M"
    )

def to_decimal_year(dt):
    year_start = datetime(dt.year, 1, 1)
    year_end = datetime(dt.year + 1, 1, 1)
    return dt.year + ((dt - year_start).total_seconds() / (year_end - year_start).total_seconds())

data['gtime'] = data['gtime'].apply(to_decimal_year)

data.to_csv("input_for_CANYONB.csv")