"""
File: calc_oxycline_features.py
-------------------------------

Estimates organic matter C:N:P ratio and AOU in the oxycline
using a robust least squares model between 45-95 dbar.

INPUT:
    :falkor_clean.csv: .csv file with cruise data and the following columns:
    'Station'
    'lon'
    'lat'
    'T'
    'S'
    'P'
    'O2'
    'sigma0'
    'rho'
    'DIC',
    'DIP'
    'NO3'
    'NO2'
    'NH4'
OUTPUTS:
    :Calculated AOU:P, DIN:P, and DIC:P in the oxycline: printout
    :Cox calculated from AOU:P, DIN:P, and DIC:P: printout
    :plots of AOU:P, DIN:P, and DIC:P in the oxycline:
"""


# Import Libraries
import numpy as np
import pandas as pd
import statsmodels.api as sm
import gsw
import matplotlib.pyplot as plt

# Load Data
falkor = pd.read_csv("falkor_clean.csv")
falkor = falkor[
    [
        "Station",
        "lon",
        "lat",
        "T",
        "S",
        "P",
        "O2",
        "sigma0",
        "rho",
        "DIC",
        "DIP",
        "NO3",
        "NO2",
        "NH4",
    ]
]
st = np.array((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19))
falkor = falkor.iloc[np.isin(falkor["Station"], st)]
print(falkor.columns)

# Calculate AOU
SA = gsw.SA_from_SP(falkor["S"], falkor["P"], falkor["lon"], falkor["lat"])
PT = gsw.pt0_from_t(SA, falkor["T"], falkor["P"])
O2sat = gsw.O2sol_SP_pt(falkor["S"], PT)
AOU = O2sat - falkor["O2"]
falkor["AOU"] = AOU
falkor["DIN"] = falkor["NO3"] + falkor["NO2"] + falkor["NH4"]

# Calculate -O2:P, AOU, C:N:P for different O2 ranges (robust fit)
# what happens if you define the oxycline using the oxygen data?
oxycline = falkor.where((falkor["P"] > 45) & (falkor["P"] < 95))
oxycline = oxycline.dropna(axis=0)

# regress O2, AOU, DIC, and DIN on phosphate (DIP) to calculate oxycline C:N:P and AOU:P
xx = sm.add_constant(oxycline["DIP"])
# Huber's t robust regression assigns varying weight values to outliers based on the
# predicted residual values, reduces the contribution of outliers to the regression,
# and increases the accuracy of slope value (Huber, 1973).
rO2 = sm.RLM(oxycline["O2"], xx, M=sm.robust.norms.HuberT()).fit()
rAOU = sm.RLM(oxycline["AOU"], xx, M=sm.robust.norms.HuberT()).fit()
rC = sm.RLM(oxycline["DIC"], xx, M=sm.robust.norms.HuberT()).fit()
rN = sm.RLM(oxycline["DIN"], xx, M=sm.robust.norms.HuberT()).fit()

# Calculate respiration quotient and Cox (oxidation state of organic carbon)
"""
this comes from calculating Cox from:
1. CaHbOc, which gives us a*Cox = 2c-b
2. Balancing the respiration reaction CaHbOc -> a CO2 + yH2O, which phrases 2c - b (and thus Cox) in terms of things we can measure:
2c - b = 4a - 4r
...where r = AOU:P and a = DIC:P, calculated from data in the oxycline
"""
Cox1 = round(
    4 - 4 * (rAOU.params[1] - 2 * rN.params[1]) / rC.params[1], 2
)  # Cox with nitrification = 4 - 4*(AOU:P - N:P)/C:P
Cox2 = round(
    4 - 4 * (rAOU.params[1]) / rC.params[1], 2
)  # Cox with no nitrification = 4 - 4*AOU:P/C:P

# Print Results
text1 = "RLS C:N:P is {}+-{} : {}+-{} : 1 (se)".format(
    round(rC.params[1], 3),
    round(rC.bse[1], 3),
    round(rN.params[1], 3),
    round(rN.bse[1], 3),
)
text2 = "RLS AOU:P is {}+-{}".format(round(rAOU.params[1], 3), round(rAOU.bse[1], 3))
# text3 should match values in Table S2
text3 = "RLS Cox calculated with AOU:P is {}, and {} with no nitrification".format(
    Cox1, Cox2
)


print(text1, text2, text3, sep="\n")


# plot regressions quickly
def quickplot(xvar, yvar, regressionparams):
    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    ax.scatter(oxycline[xvar], oxycline[yvar])
    xfit = np.linspace(oxycline[xvar].min(), oxycline[xvar].max())
    ax.plot(xfit, xfit * regressionparams[1] + regressionparams[0], color="k")
    ax.set_xlabel(xvar)
    ax.set_ylabel(yvar)
    plt.tight_layout()
    plt.show()


quickplot("DIP", "O2", rO2.params)
quickplot("DIP", "AOU", rAOU.params)
quickplot("DIP", "DIC", rC.params)
quickplot("DIP", "DIN", rN.params)
