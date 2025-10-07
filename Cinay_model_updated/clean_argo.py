"""
File: clean_argo.py
---------------------

Selects the measurement with 2 and 6 flags.
Calculates DIC using PyCO2SYS and density using TEOS-10 GSW packages.
Saves the analyzed variables and data as a .csv file.

INPUT:
    :bottle_data.csv: .csv file with cruise data and the following columns:
        'STNNBR'
        'CASTNO'
        'LATITUDE'
        'LONGITUDE'
        'CTDPRS'
        'CTDTMP'
        'CTDSAL'
        'CTDS_FLAG_W'
        'CTDOXY'
        'CTDOXY_FLAG_W'
        'FLOR'
        'FLOR_FLAG_W'
        'NOx_FLAG_W'
        'NITRIT_BabLab'
        'NITRIT_FLAG_W'
        'PHSPHT'
        'PHSPHT_FLAG_W'
        'NH4_FLAG_W'
        'PH_TOT'
        'pH_insitu_flag_W'
        'PH_TEMP'
        'TA_FLAG_W'
OUTPUT:
    :data_clean.csv: .csv file with bad flagged data filtered out and calculated params:
        -absolute salinity (calculated with gsw)
        -conservative temperature (calculated with gsw)
        -DIC (calculated with pyCO2SYS)
        -pH insitu (calculated with pyCO2SYS)
"""

## Import Libraries
import pandas as pd
import numpy as np
import gsw

data = pd.read_csv("CANYONB_output.csv")

# remove extra index column
del data["Var1"]

# rename columns
cols = {
    'Station': "station",
    'Lat [°N]': "lat",
    'Lon [°E]': "lon",
    'pres': "press",
    'temp': "temperature",
    'QF_2':'temp_flag',
    'psal': "sal",
    'QF_3': "sal_flag",
    'doxy': "O2",
    'QF_6': "O2_flag",
    'Nitrate__mol_kg_':"NO3",
    'QF_8':'NO3_flag',
    'Nitrite__mol_kg_': "NO2",
    'QF_17': "NO2_flag",
    'PO4_CANY': "phosphate",
    'pHinsitu_Total_': "pH_insitu",
    'QF_12': "pH_insitu_flag",
    'TALK_LIAR__mol_kg_':'TA',
    'QF_14':"TA_flag",
    'DIC_LIAR__mol_kg_':'DIC',
    'QF_15':'DIC_flag'
}
data = data.rename(columns=cols)

# Select good flagged data

flag_good = np.array([0]) # 0 = good, 8 = bad
idx_flag = np.where(
    (np.isin(data["temp_flag"], flag_good))
    & (np.isin(data["sal_flag"], flag_good))
    & (np.isin(data["O2_flag"], flag_good))
    & (np.isin(data["NO3_flag"], flag_good))
    & (np.isin(data["pH_insitu_flag"], flag_good))
    & (np.isin(data["TA_flag"], flag_good))
    & (np.isin(data["NO2_flag"], flag_good))
)

# calculated parameters
SA = gsw.SA_from_SP(data["sal"], data["press"], data["lon"], data["lat"])
CT = gsw.CT_from_t(SA, data["temperature"], data["press"])
data["sigma0"] = gsw.density.sigma0(SA, CT)
data["rho"] = gsw.density.rho(SA, CT, data["press"])

data["Nstar"] = ((data.NO2 + data.NO3) - 16 * data.phosphate + 2.9)

# optional - Calculate Carbonate System with pyCO2sys instead of using LIAR vars
'''
Z = pyco2.sys(
    par1=data["TA"],
    par2=data["pH_tot"],
    par1_type=1,
    par2_type=3,
    salinity=data["sal"],
    temperature=25,
    temperature_out=data["temperature"],
    pressure=0,
    pressure_out=data["press"],
    #total_phosphate=data["phosphate"],
    #total_ammonia=data["NH4"],
    opt_pH_scale=1,
    opt_k_carbonic=10,
    opt_total_borate=2,
    opt_k_fluoride=2,
)
data["DIC"] = Z["dic"]
data["pH_insitu"] = Z["pH_out"]
'''

# Remove Bad Data and Save New Vectors
depth = np.array(data['Depth_m_'])[idx_flag]  # degrees C
date = np.array(data['mon_day_yr'])[idx_flag]
T = np.array(data["temperature"])[idx_flag]  # degrees C
S = np.array(data["sal"])[idx_flag]  # psu
P = np.array(data["press"])[idx_flag]  # dbar
rho = np.array(data["rho"])[idx_flag]  # kg/m3
sigma0 = np.array(data["sigma0"])[idx_flag]  # kg/m3 #
DIC = np.array(data["DIC"])[idx_flag]  # umol/kg
DIP = np.array(data["phosphate"])[idx_flag] # umol/kg
NO2 = np.array(data["NO2"])[idx_flag] # umol/kg
NO3 = np.array(data["NO3"])[idx_flag] # umol/kg
Nstar = ((NO2 + NO3) - 16 * DIP + 2.9)  # umol/kg  ## Change the DIP coefficient (either 11.4 or 16) for sensitivity analyses.
TA = np.array(data["TA"])[idx_flag]  # umol/kg
pH = np.array(data["pH_insitu"])[idx_flag]
O2 = np.array(data["O2"])[idx_flag]  # umol/kg
station = np.array(data["station"])[idx_flag]
lat = np.array(data["lat"])[idx_flag]
lon = np.array(data["lon"])[idx_flag]
#pH_tot = np.array(data["pH_tot"])[idx_flag]
#pH_temp = np.array(data["pH_temp"])[idx_flag]
#fluor = np.array(data["fluor"])[idx_flag]
data_array = np.array(
    [
        depth,
        date,
        station,
        #cast,
        lat,
        lon,
        T,
        S,
        P,
        rho,
        sigma0,
        DIC,
        DIP,
        NO3,
        NO2,
        #NH4,
        Nstar,
        TA,
        pH,
        #pH_tot,
        #pH_temp,
        O2,
        #fluor,
    ]
).T

data_clean = pd.DataFrame(
    data=data_array,
    columns=[
        "Depth",
        "Date",
        "Station",
        #"cast",
        "lat",
        "lon",
        "T",
        "S",
        "P",
        "rho",
        "sigma0",
        "DIC",
        "DIP",
        "NO3",
        "NO2",
        #"NH4",
        "Nstar",
        "TA",
        "pH",
        #"pH_tot",
        #"pH_temp",
        "O2",
        #"fluor",
    ],
)

# Save the clean version
data_clean.to_csv("data_clean.csv")
