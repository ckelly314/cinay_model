"""
File: calc_argo_outputs.py
----------------------------

Calculates the relative reaction rates and relative contributions based on the R matrix
input selected with Monte Carlo simuations. For the delta_tracer values, the script fits
the tracer data from data_clean.csv file within each layer (defined in the script).
Saves relative contributions as a text file for further plotting.

Inputs:
    :data_clean.csv: .csv file with cruise data and the following columns:
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
    :R.txt: comma-delimited text file containing R matrix
Outputs:
    :layers.csv: 17x1 matrix, defining the sublayers along which we apply calculations.
        Each row represents a density surface.
    :reaction_matrix.csv: 6x5 matrix ("R matrix") containing a theoretical Δtracer:DIC for different reactions.
        Each row represents a tracer (n=6), and each column represents a reaction (n=5).
    :slopes_mean.csv: 16x6 matrix, where:
        Each row represents a sublayer, bounded on each side by the density surfaces defined in "sublayers".
        Each column represents a tracer.
        Each cell is the measured slope of a given tracer against DIC in a given density layer.
    :slopes_se.csv: same as slopes_mean.csv, but standard errors of slopes
    :coeffs_mean.csv: 16x5 matrix where each row represents a sublayer and each column represents a process.
        The values in each row are the mean values of the X matrix from a Monte Carlo simulation with 1000 iterations.
    :coeffs_se.csv: same as coeffs_mean but standard errors of Monte Carlo simulation.
    :relative_importances_mean.csv: 16x8 matrix, where:
        Each row represents a sublayer, and is the mean of 1000 Monte Carlo iterations
        Each column represents a reaction: anmx, denit, ox, red, nitox, dnrn, caco3, otherDIC
    :relative_importances_se.csv: same as relative_importances_mean but with MC standard errors
    :relative_importances_nitrite_mean.csv: 16x4 matrix, where:
        Each row represents a sublayer, and is the mean of 1000 Monte Carlo iterations
        Each column represents a reaction: anmx, denit, nitox, dnrn
    :relative_importances_nitrite_mean.csv: same as relative_importances_nitrite_mean but with MC standard errors
    :residuals.csv: 16x6 matrix containing (measured Δtracer:DIC) - (theoretical Δtracer:DIC)
        Each row is a density layer;
        Each column is a tracer: ΔNO3-, ΔNO2-, ΔNH4+, ΔN*, ΔTA, and ΔDIC
    :percentage_residuals.csv: residuals as a percentage of the observed slopes
"""

## Import Libraries
import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.optimize as sc
import matplotlib.pyplot as plt
import seaborn as sns
import warnings # suppress warning for mean of empty slice
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

### USER INPUTS ###
# stations = np.arange(1,15,1, dtype=float) # Stations selected for analysis
#stations = np.arange(1,57) #np.array((1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14))
#layers = np.array(
    #[25, 25.5, 25.8, 26.21, 26.31, 26.44, 26.67, 27, 27.3], dtype=float
#    [25.5, 26, 27.3], dtype=float
#)  # in sigma0
chunkID = "lowpoclownitrite"
makeplots = True
if chunkID == 1:
    stations = np.arange(1,21)
    # Define Layers
    layers = np.array([25.8, 26.3, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == 2:
    stations = np.arange(21,40)
    # Define Layers
    layers = np.array([25.75, 26.35, 26.55,27.2], dtype=float)  # in sigma0
elif chunkID == 3:
    stations = np.arange(40,44)
    # Define Layers
    layers = np.array([25.75, 26.35,26.5,27.3], dtype=float)  # in sigma0
elif chunkID == 4:
    stations = np.arange(44,48)
    # Define Layers
    layers = np.array([25.75, 26.35,26.55,27.3], dtype=float)  # in sigma0
elif chunkID == 5:
    stations = np.arange(48,57)
    # Define Layers
    layers = np.array([25.75, 26.2,26.55,27.3], dtype=float)  # in sigma0\
elif chunkID == 6:
    stations = np.arange(58,71)
    # Define Layers
    layers = np.array([25.99, 26.25, 26.55,27.2], dtype=float)  # in sigma0
elif chunkID == 7:
    stations = np.arange(1,72)
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == "highnitrite":
    stations = np.arange(1,30)
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == "lownitrite":
    stations = np.arange(30,72)
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == "highpoc":
    stations = np.arange(20,28)
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == "lowpochighnitrite":
    stations = np.concatenate((np.arange(1,20), np.arange(28,30)))
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
elif chunkID == "lowpoclownitrite":
    stations = np.arange(30,72)
    # Define Layers
    layers = np.array([25.8, 26.1, 26.35, 26.5,27.2], dtype=float)  # in sigma0
###################

def quickplot(xvar, yvar, cvar, huberregression, olsregression, xlabel, ylabel, lower_boundary, upper_boundary,
             huber_loss, sse):
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    cax = ax.scatter(xvar, yvar,
        c=cvar)
    xfit = np.linspace(xvar.min(), xvar.max())
    
    ax.plot(xfit, xfit * huberregression.params[1] + huberregression.params[0], color="k",
        label = "Robust HuberT")
    ax.plot(xfit, xfit * olsregression.params[1] + olsregression.params[0], color="k", linestyle = ":",
        label = "OLS")
    
    textstr = f"Huber Loss={huber_loss:.4}\nSSE={sse:.4}\np={huberregression.pvalues[1]:.3}"
    
    ax.text(1.5, 0.0, textstr, transform = ax.transAxes,
        verticalalignment = "bottom")
    
    fig.colorbar(cax).set_label(r"$\sigma_{\theta}$")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(fr'$\sigma_{{\theta}}${lower_boundary:.4}-{upper_boundary:.4}')
    ax.legend(bbox_to_anchor=(1.5,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"figures/chunk{chunkID}/{lower_boundary}-{upper_boundary}_{ylabel}vs{xlabel}.pdf",)
    plt.show()

def plotobsest(xvar, yvar, cvar, regression, xlabel, ylabel, lower_boundary, upper_boundary, slope_est):
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    cax = ax.scatter(xvar, yvar,
        c=cvar)
    xfit = np.linspace(xvar.min(), xvar.max())
    ax.plot(xfit, xfit * regression.params[1] + regression.params[0], color="k",
        label = "obs")
    y1 = np.median(xfit) * regression.params[1] + regression.params[0]
    y2 =  np.median(xfit) * slope_est
    offset = y1-y2
    ax.plot(xfit, xfit * slope_est + offset, color="r",
        label = "est")
    textstr = f"p={regression.pvalues[1]:.2}"
    ax.text(0.05, 0.95, textstr, transform = ax.transAxes,
        verticalalignment = "top")
    fig.colorbar(cax).set_label(r"$\sigma_{\theta}$")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(fr'$\sigma_{{\theta}}${lower_boundary:.4}-{upper_boundary:.4}')
    ax.legend(bbox_to_anchor=(1.4,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"figures/chunk{chunkID}/{lower_boundary}-{upper_boundary}_{ylabel}vs{xlabel}.pdf",)
    plt.show()

def huberregression(Y, xx):
    rf = sm.RLM(Y, xx, M=sm.robust.norms.HuberT()).fit()
    # Residuals
    resid = rf.resid
    delta = sm.robust.norms.HuberT().t
    
    # Compute Huber loss manually
    huber_loss = np.where(np.abs(resid) <= delta,
                          0.5 * resid**2,
                          delta * (np.abs(resid) - 0.5 * delta))
    
    huber_loss = np.sum(huber_loss)
    
    print(f"Total Huber loss: {huber_loss:.3f}")

    ols_model = sm.OLS(Y, xx).fit()
    sse = np.sum((Y - ols_model.fittedvalues) ** 2)
    print(f"OLS SSE: {sse:.3f}")
    
    return (rf, ols_model, huber_loss, sse)

# Set Directory
fpath = "output/OM_variations/experimental2/{}"
#fpath = 'output/OM_variations/anderson/{}'

# Import Data
data = pd.read_csv("data_clean.csv")
data = data.dropna()
R = np.loadtxt(fpath.format("R.txt"), delimiter=",")
Rsolve = R[[0,1,3,4,5]] # remove row for NH4

# check dates
stationdict = {"highnitrite":np.arange(1,30),
               "lownitrite":np.arange(30,72),
               "highpoc":np.arange(20,28),
               "lowpochighnitrite":np.concatenate((np.arange(1,20), np.arange(28,30))),
               "lowpoclownitrite":np.arange(30,72)}  # in sigma0
mindate = data[np.isin(data.Station, stationdict[chunkID])].Date.unique().min()
maxdate = data[np.isin(data.Station, stationdict[chunkID])].Date.unique().max()
print(chunkID, mindate, maxdate)

# Define Inputs
divider = 2  # Number of sublayers in each layer
K = 10000  # Number of Iterations for Monte Carlo Error Propagation

### Data Preparation ###
# Select Stations
idx_station = np.where(np.isin(data["Station"], stations))

# Save Data Vectors
rho = np.array(data["rho"])[idx_station]  # kg/m3
sigma0 = np.array(data["sigma0"])[idx_station]  # kg/m3
DIC = np.array(data["DIC"])[idx_station]  # umol/kg
DIP = np.array(data["DIP"])[idx_station]  # umol/kg
NO2 = np.array(data["NO2"])[idx_station]  # umol/kg
NO3 = np.array(data["NO3"])[idx_station]  # umol/kg
#NH4 = np.array(data["NH4"])[idx_station]  # umol/kg
Nstar = np.array(data["Nstar"])[idx_station]  # umol/kg
TA = np.array(data["TA"])[idx_station]  # umol/kg
pH = np.array(data["pH"])[idx_station]
O2 = np.array(data["O2"])[idx_station]  # umol/kg
depth = np.array(data["Depth"])[idx_station]

# divide layers up into sublayers
sl = np.zeros((len(layers) - 1, divider + 1))
for i in np.arange(0, len(layers) - 1):
    sl[i,] = np.linspace(layers[i], layers[i + 1], divider + 1)
sublayers = np.unique(sl)

### Robust Linear Regression ###
# Create Results Arrays
# for each slice of water between defined sublayers, store (slope, slope error, intercept, intercept error)
fitting_NO3 = np.zeros((len(sublayers) - 1, 4))
fitting_NO2 = np.zeros((len(sublayers) - 1, 4))
#fitting_NH4 = np.zeros((len(sublayers) - 1, 4))
fitting_Nstar = np.zeros((len(sublayers) - 1, 4))
fitting_TA = np.zeros((len(sublayers) - 1, 4))
fitting_DIC = np.zeros((len(sublayers) - 1, 4))

# Run Robost Linear Regression
# store number of data points in each layer for model stats
npoints = np.zeros((len(sublayers) - 1,1))
sigmas = np.zeros((len(sublayers) - 1,1))
depths = np.zeros((len(sublayers) - 1,1))
depthserr = np.zeros((len(sublayers) - 1,1))

for i in np.arange(0, len(sublayers) - 1):
    lower_boundary = sublayers[i]
    upper_boundary = sublayers[i + 1]
    idx_layer = np.where((sigma0 >= lower_boundary) & (sigma0 <= upper_boundary))
    if chunkID == 5:
        idx_layer = np.where((sigma0 >= lower_boundary) & (sigma0 <= upper_boundary)
            & (DIC < 2300))

    # Robust Regression
    if len(DIC[idx_layer]) > 0:
        # DIC is the x-variable against which we regress NO2-, NO2-, N*, TA, and DIC (the last of which should be 1:1)
        xx = sm.add_constant(DIC[idx_layer])
        npoints[i] = len(xx)
        sigmas[i] = (lower_boundary + upper_boundary)/2
        depths[i] = np.nanmean(depth[idx_layer])
        depthserr[i] = np.nanstd(depth[idx_layer])

        # robust regression
        rf_NO3, ols_NO3, loss_NO3, sse_NO3 = huberregression(NO3[idx_layer], xx)
        rf_NO2, ols_NO2, loss_NO2, sse_NO2 = huberregression(NO2[idx_layer], xx)
        #rf_NH4, ols_NH4, loss_NH4, sse_NH4 = huberregression(NH4[idx_layer], xx)
        rf_Nstar, ols_Nstar, loss_Nstar, sse_Nstar = huberregression(Nstar[idx_layer], xx)
        rf_TA, ols_TA, loss_TA, sse_TA = huberregression(TA[idx_layer], xx)
        rf_DIC, ols_DIC, loss_DIC, sse_DIC = huberregression(DIC[idx_layer], xx)

        # Save Results
        # columns in "fitting_NO3" are slope, slope err, int, int err
        # rows correspond to each slice of water between sublayers
        fitting_NO3[i, :] = np.array(
            [rf_NO3.params[1], rf_NO3.bse[1], rf_NO3.params[0], rf_NO3.bse[0]]
        )
        fitting_NO2[i, :] = np.array(
            [rf_NO2.params[1], rf_NO2.bse[1], rf_NO2.params[0], rf_NO2.bse[0]]
        )
        #fitting_NH4[i, :] = np.array(
        #    [rf_NH4.params[1], rf_NH4.bse[1], rf_NH4.params[0], rf_NH4.bse[0]]
        #)
        fitting_Nstar[i, :] = np.array(
            [rf_Nstar.params[1], rf_Nstar.bse[1], rf_Nstar.params[0], rf_Nstar.bse[0]]
        )
        fitting_TA[i, :] = np.array(
            [rf_TA.params[1], rf_TA.bse[1], rf_TA.params[0], rf_TA.bse[0]]
        )
        fitting_DIC[i, :] = np.array(
            [rf_DIC.params[1], rf_DIC.bse[1], rf_DIC.params[0], rf_DIC.bse[0]]
        )

        quickplot(DIC[idx_layer], NO3[idx_layer], sigma0[idx_layer], rf_NO3, ols_NO3,
            "[DIC]", "[NO3-]", lower_boundary, upper_boundary, loss_NO3, sse_NO3)
        quickplot(DIC[idx_layer], NO2[idx_layer], sigma0[idx_layer], rf_NO2, ols_NO2,
            "[DIC]", "[NO2-]", lower_boundary, upper_boundary, loss_NO2, sse_NO2)
        quickplot(DIC[idx_layer], Nstar[idx_layer], sigma0[idx_layer], rf_Nstar, ols_Nstar,
            "[DIC]", "N*", lower_boundary, upper_boundary, loss_Nstar, sse_Nstar)
        quickplot(DIC[idx_layer], TA[idx_layer], sigma0[idx_layer], rf_TA, ols_TA,
            "[DIC]", "TA", lower_boundary, upper_boundary, loss_TA, sse_TA)
        quickplot(DIC[idx_layer], DIC[idx_layer], sigma0[idx_layer], rf_DIC, ols_DIC,
            "[DIC]", "[DIC]", lower_boundary, upper_boundary, loss_DIC, sse_DIC)


    else:
        print(f"no data {sublayers[i]}-{sublayers[i + 1]}")


# Save Slopes
slopes_mean = np.array( # organize slopes together into one array
    [
        fitting_NO3[:, 0],
        fitting_NO2[:, 0],
        #fitting_NH4[:, 0],
        fitting_Nstar[:, 0],
        fitting_TA[:, 0],
        fitting_DIC[:, 0],
    ]
).T
slopes_se = np.array( # organize slope errors into another array
    [
        fitting_NO3[:, 1],
        fitting_NO2[:, 1],
        #fitting_NH4[:, 1],
        fitting_Nstar[:, 1],
        fitting_TA[:, 1],
        fitting_DIC[:, 1],
    ]
).T

# Layer Labels
ylabels = []
for i in range(len(sublayers)-1):
    label = f"sigma_theta = {sublayers[i]:.4} - {sublayers[i+1]:.4}"
    ylabels.append(label)

# generate Excel file for visual inspection
slopeoutput = pd.DataFrame(slopes_mean, columns = ["NO3vsDICslope", "NO2vsDICslope", "NstarvsDICslope","TAvsDICslope","DICvsDICslope"])
slopeoutput2 = pd.DataFrame(slopes_se, columns = ["NO3vsDICerr", "NO2vsDICerr", "NstarvsDICerr","TAvsDICerr","DICvsDICerr"])
slopeoutput["layer"] = ylabels
slopeoutput2["layer"] = ylabels
slopeoutput = slopeoutput.set_index("layer", drop=True)
slopeoutput2 = slopeoutput2.set_index("layer", drop=True)
slopeoutput = slopeoutput.join(slopeoutput2)
slopeoutput["n"] = npoints
slopeoutput["sigma0"] = sigmas
slopeoutput["depth"] = depths
slopeoutput["deptherr"] = depthserr
slopeoutput = slopeoutput[["n","sigma0","depth","deptherr",
    'NO3vsDICslope','NO3vsDICerr',
                           'NO2vsDICslope','NO2vsDICerr',
                           'NstarvsDICslope','NstarvsDICerr',
                           'TAvsDICslope','TAvsDICerr',
                           'DICvsDICslope','DICvsDICerr']]
slopeoutput.to_excel(f"output/chunk{chunkID}/slopesDF.xlsx")

### Reaction Coefficient, Relative Contributions, and Residuals Calculations with Monte Carlo Error Propagation ###
# Create Results Arrays
slope_iter = np.zeros((slopes_mean.shape[1]))
coeff_iter = np.zeros((K, 5))
coeff_mean = np.zeros((len(sublayers) - 1, 5))
coeff_se = np.zeros((len(sublayers) - 1, 5))
relimp_mean = np.zeros((len(sublayers) - 1, 8))
relimp_se = np.zeros((len(sublayers) - 1, 8))
relnit_mean = np.zeros((len(sublayers) - 1, 4))
relnit_se = np.zeros((len(sublayers) - 1, 4))
relimp_iter = np.zeros((K, 8))
anmx_iter = np.zeros((K))
denit_iter = np.zeros((K))
ox_iter = np.zeros((K))
red_iter = np.zeros((K))
nitox_iter = np.zeros((K))
dnrn_iter = np.zeros((K))
caco3_iter = np.zeros((K))
otherDIC_iter = np.zeros((K))
denit2_iter = np.zeros((K))
relnit_anmx_iter = np.zeros((K))
relnit_denit_iter = np.zeros((K))
relnit_nitox_iter = np.zeros((K))
relnit_dnrn_iter = np.zeros((K))
relnit_iter = np.zeros((K, 4))
slopes_obs = np.zeros((len(sublayers) - 1, slopes_mean.shape[1]))
slopes_est = np.zeros((len(sublayers) - 1, slopes_mean.shape[1]))
residuals = np.zeros((len(sublayers) - 1, slopes_mean.shape[1]))
residuals_perc = np.zeros((len(sublayers) - 1, slopes_mean.shape[1]))
nstar_slope = np.zeros((K))

# Calculate Coeffs, Relative Contributions, and Residuals
for i in np.arange(0, len(sublayers) - 1):  # do Monte Carlo simulation for each layer
    # Run Monte Carlo
    for k in range(K):  # K is the number of iterations
        for j in range(
            slopes_mean.shape[1]
        ):  # pick a value for measured Δtracer:DIC based on normal distribution
            slope_iter[j] = np.random.normal(
                loc=slopes_mean[i, j], scale=slopes_se[i, j], size=1
            )
        # Calculate Coeffs
        # nnls = "non-negatcie"
        # we use nnls instead of np.linalg.solve because we don't want any negative numbers
        c_temp, rnorm = sc.nnls(Rsolve, slope_iter)  # solve for Χ matrix
        coeff_iter[k, :] = c_temp / np.sum(c_temp)  # make sure X matrix sums to 1

        nstar_slope[k] = np.random.normal(
           loc=slopes_mean[i, 4], scale=slopes_se[i, 4], size=1
        )
        # Calculate Relative Importances
        ''' OLD CALCS
        if abs(R[3, 2] * c_temp[2] + R[3, 1] * c_temp[1]) != 0:
            anmx_iter[k] = ((R[3, 2] * c_temp[2]) / (R[3, 2] * c_temp[2] + R[3, 1] * c_temp[1]) * 100)
            denit_iter[k] = ((R[3, 1] * c_temp[1]) / (R[3, 2] * c_temp[2] + R[3, 1] * c_temp[1]) * 100)
        else:
            anmx_iter[k] = 0
            denit_iter[k] = 0
        '''

        #if abs(R[1, 2] * coeff_iter[k, :][2] + R[1, 1] * coeff_iter[k, :][1]) != 0:
        if (coeff_iter[k, :][2] != 0) | (coeff_iter[k, :][1] != 0): # if relative reaction rate of anammox OR denit != 0,
            # calculate proportional contributions of anammox and denit to nitrite drawdown
            anmx_iter[k] = abs(R[1, 2] * coeff_iter[k, :][2]) / (abs(R[1, 2] * coeff_iter[k, :][2]) + abs(R[1, 1] * coeff_iter[k, :][1])) * 100
            denit_iter[k] = abs(R[1, 1] * coeff_iter[k, :][1]) / (abs(R[1, 2] * coeff_iter[k, :][2]) + abs(R[1, 1] * coeff_iter[k, :][1]))* 100
        else: # otherwise, set relative proportions of anammox and denit to NaN
            anmx_iter[k] = np.nan
            denit_iter[k] = np.nan

        # calculate net nitrite oxidation 
        ox_iter[k] = ( # little bit of nitrite is effectively oxidized to nitrate during anammox
            (R[1, 3] * coeff_iter[k, 3] - R[0, 2] * coeff_iter[k, 2]) # R[NO2, nitrox]*coeffs[nitrox] - R[NO3, anammox]*coeffs[anammox]
            / (
                R[1, 2] * coeff_iter[k, 2]
                + R[1, 1] * coeff_iter[k, 1]
                + R[1, 3] * coeff_iter[k, 3]
                - R[0, 2] * coeff_iter[k, 2]
            )
            * 100
        )
        red_iter[k] = (
            (R[1, 2] * coeff_iter[k, 2] + R[1, 1] * coeff_iter[k, 1])
            / (
                R[1, 2] * coeff_iter[k, 2]
                + R[1, 1] * coeff_iter[k, 1]
                + R[1, 3] * coeff_iter[k, 3]
                - R[0, 2] * coeff_iter[k, 2]
            )
            * 100
        )
        nitox_iter[k] = (
            abs(R[1, 3] * coeff_iter[k, 3])
            / abs(R[1, 3] * coeff_iter[k, 3] - R[1, 0] * coeff_iter[k, 0])
            * 100
        )
        dnrn_iter[k] = (
            abs(R[1, 0] * coeff_iter[k, 0])
            / abs(R[1, 3] * coeff_iter[k, 3] - R[1, 0] * coeff_iter[k, 0])
            * 100
        )
        caco3_iter[k] = (
            abs(coeff_iter[k, 4])
            / (abs(coeff_iter[k, 0] + coeff_iter[k, 1] - coeff_iter[k, 2] - coeff_iter[k, 3]) + coeff_iter[k, 4])
            * 100
        )
        otherDIC_iter[k] = (
            abs(coeff_iter[k, 0] + coeff_iter[k, 1] - coeff_iter[k, 2] - coeff_iter[k, 3])
            / (abs(coeff_iter[k, 0] + coeff_iter[k, 1] - coeff_iter[k, 2] - coeff_iter[k, 3]) + coeff_iter[k, 4])
            * 100
        )
        relimp_iter[k, :] = [
            anmx_iter[k],
            denit_iter[k],
            ox_iter[k],
            red_iter[k],
            nitox_iter[k],
            dnrn_iter[k],
            caco3_iter[k],
            otherDIC_iter[k],
        ]
        relnit_anmx_iter[k] = (
            abs(R[1, 2] * coeff_iter[k, 2])
            / abs(
                -R[1, 2] * coeff_iter[k, 2]
                - R[1, 1] * coeff_iter[k, 1]
                - R[1, 3] * coeff_iter[k, 3]
                + R[1, 0] * coeff_iter[k, 0]
            )
            * 100
        )
        relnit_denit_iter[k] = (
            abs(R[1, 1] * coeff_iter[k, 1])
            / abs(
                -R[1, 2] * coeff_iter[k, 2]
                - R[1, 1] * coeff_iter[k, 1]
                - R[1, 3] * coeff_iter[k, 3]
                + R[1, 0] * coeff_iter[k, 0]
            )
            * 100
        )
        relnit_nitox_iter[k] = (
            abs(R[1, 3] * coeff_iter[k, 3])
            / abs(
                -R[1, 2] * coeff_iter[k, 2]
                - R[1, 1] * coeff_iter[k, 1]
                - R[1, 3] * coeff_iter[k, 3]
                + R[1, 0] * coeff_iter[k, 0]
            )
            * 100
        )
        relnit_dnrn_iter[k] = (
            abs(R[1, 0] * coeff_iter[k, 0])
            / abs(
                -R[1, 2] * coeff_iter[k, 2]
                - R[1, 1] * coeff_iter[k, 1]
                - R[1, 3] * coeff_iter[k, 3]
                + R[1, 0] * coeff_iter[k, 0]
            )
            * 100
        )
        relnit_iter[k, :] = [
            relnit_anmx_iter[k],
            relnit_denit_iter[k],
            relnit_nitox_iter[k],
            relnit_dnrn_iter[k],
        ]

    # save out full Monte Carlo array so we can do statistical tests with output
    np.savetxt(
        #fpath.format("slopes_mean.csv"),
        f"output/chunk{chunkID}/layer{i}_coeff_iter.csv",
        coeff_iter,
        delimiter=",",
        header="DNRN, Denitrification, Anammox, Nitrite Oxidation, CaCO3 Dissolution ", #"NO3, NO2, NH4, Nstar, TA, DIC",
    )
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning) # suppress warnings for mean of empty slice
        coeff_mean[i, :] = np.mean(
            coeff_iter, axis=0
        )  # take the mean of the 1000 Monte Carlo iterations
        coeff_se[i, :] = np.std(coeff_iter, axis=0)  # /math.sqrt(K)
        relimp_mean[i, :] = np.nanmean(relimp_iter, axis=0)
        relimp_se[i, :] = np.nanstd(relimp_iter, axis=0)  # /math.sqrt(K)
        relnit_mean[i, :] = np.nanmean(relnit_iter, axis=0)
        relnit_se[i, :] = np.nanstd(relnit_iter, axis=0)  # /math.sqrt(K)

        # Calculate Residuals
        slopes_obs[i, :] = slopes_mean[i, :].T
        slopes_est[i, :] = np.dot(Rsolve, coeff_mean[i, :])
        residuals[i, :] = slopes_obs[i, :] - slopes_est[i, :]
        residuals_perc[i, :] = abs(residuals[i, :]) / abs(slopes_obs[i, :]) * 100

if makeplots == True:
    for i in np.arange(0, len(sublayers) - 1):
        lower_boundary = sublayers[i]
        upper_boundary = sublayers[i + 1]
        idx_layer = np.where((sigma0 >= lower_boundary) & (sigma0 <= upper_boundary))
        if chunkID == 5:
            idx_layer = np.where((sigma0 >= lower_boundary) & (sigma0 <= upper_boundary)
                & (DIC < 2300))

        # Robust Regression
        if len(DIC[idx_layer]) > 0:
            xx = sm.add_constant(DIC[idx_layer])
            rf_NO3, ols_NO3, loss_NO3, sse_NO3 = huberregression(NO3[idx_layer], xx)
            rf_NO2, ols_NO2, loss_NO2, sse_NO2 = huberregression(NO2[idx_layer], xx)
            #rf_NH4, ols_NH4, loss_NH4, sse_NH4 = huberregression(NH4[idx_layer], xx)
            rf_Nstar, ols_Nstar, loss_Nstar, sse_Nstar = huberregression(Nstar[idx_layer], xx)
            rf_TA, ols_TA, loss_TA, sse_TA = huberregression(TA[idx_layer], xx)
            rf_DIC, ols_DIC, loss_DIC, sse_DIC = huberregression(DIC[idx_layer], xx)

            plotobsest(DIC[idx_layer], NO3[idx_layer], sigma0[idx_layer], rf_NO3, 
                "[DIC]", "[NO3-]", lower_boundary, upper_boundary, slopes_est[i,0])
            plotobsest(DIC[idx_layer], NO2[idx_layer], sigma0[idx_layer], rf_NO2, 
                "[DIC]", "[NO2-]", lower_boundary, upper_boundary, slopes_est[i,1])
            plotobsest(DIC[idx_layer], Nstar[idx_layer], sigma0[idx_layer], rf_Nstar, 
                "[DIC]", "N*", lower_boundary, upper_boundary, slopes_est[i,2])
            plotobsest(DIC[idx_layer], TA[idx_layer], sigma0[idx_layer], rf_TA, 
                "[DIC]", "TA", lower_boundary, upper_boundary, slopes_est[i,3])
            plotobsest(DIC[idx_layer], DIC[idx_layer], sigma0[idx_layer], rf_DIC, 
                "[DIC]", "[DIC]", lower_boundary, upper_boundary, slopes_est[i,4])

### Save the Outputs as a CSV file ###
np.savetxt(
    #fpath.format("slopes_mean.csv"),
    f"output/chunk{chunkID}/slopes_mean.csv",
    slopes_mean,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("slopes_se.csv"),
    f"output/chunk{chunkID}/slopes_se.csv",
    slopes_se,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("coeffs_mean.csv"),
    f"output/chunk{chunkID}/coeffs_mean.csv",
    coeff_mean,
    delimiter=",",
    header="DNRN, Anammox, Denitrification, Nitrite Oxidation, CaCO3 Dissolution ",
)
np.savetxt(
    #fpath.format("coeffs_se.csv"),
    f"output/chunk{chunkID}/coeffs_se.csv",
    coeff_se,
    delimiter=",",
    header="DNRN, Anammox, Denitrification, Nitrite Oxidation, CaCO3 Dissolution ",
)
np.savetxt(
    #fpath.format("relative_importances_mean.csv"),
    f"output/chunk{chunkID}/relative_importances_mean.csv",
    relimp_mean,
    delimiter=",",
    header="anmx, denit, ox, red, nitox, dnrn, caco3 diss, other DIC",
)
np.savetxt(
    #fpath.format("relative_importances_nitrite_mean.csv"),
    f"output/chunk{chunkID}/relative_importances_nitrite_mean.csv",
    relnit_mean,
    delimiter=",",
    header="anmx, denit, nitox, dnrn",
)
np.savetxt(
    #fpath.format("relative_importances_nitrite_se.csv"),
    f"output/chunk{chunkID}/relative_importances_nitrite_se.csv",
    relnit_se,
    delimiter=",",
    header="anmx, denit, nitox, dnrn",
)
np.savetxt(
    #fpath.format("relative_importances_se.csv"),
    f"output/chunk{chunkID}/relative_importances_se.csv",
    relimp_se,
    delimiter=",",
    header="anmx, denit, ox, red, nitox, dnrn, caco3 diss, other DIC",
)
np.savetxt(
    #fpath.format("residuals.csv"),
    f"output/chunk{chunkID}/slopes_obs.csv",
    slopes_obs,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("residuals.csv"),
    f"output/chunk{chunkID}/slopes_est.csv",
    slopes_est,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("residuals.csv"),
    f"output/chunk{chunkID}/residuals.csv",
    residuals,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("percentage_residuals.csv"),
    f"output/chunk{chunkID}/percentage_residuals.csv",
    residuals_perc,
    delimiter=",",
    header="NO3, NO2, Nstar, TA, DIC", #"NO3, NO2, NH4, Nstar, TA, DIC",
)
np.savetxt(
    #fpath.format("reaction_matrix.csv"),
    f"output/chunk{chunkID}/reaction_matrix.csv",
    R,
    delimiter=",",
    header="Columns: Reactions, Rows: Tracers",
)
np.savetxt(fpath.format("layers.csv"), sublayers, delimiter=",", header="Layers")
print("done")

