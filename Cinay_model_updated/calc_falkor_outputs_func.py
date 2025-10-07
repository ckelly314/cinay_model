"""
File: calc_falkor_outputs_func.py
---------------------------------

Turn calc_falkor_outputs.py into a function to iteratively calculate outputs
for all organic matter scenarios. - CLK

Calculates the relative reaction rates and relative contributions based on the R matrix
input selected with Monte Carlo simuations. For the delta_tracer values, the script fits
the tracer data from falkor_clean.csv file within each layer (defined in the script).
Saves relative contributions as a text file for further plotting.

Inputs:
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


def calc_falkor_outputs(fpath):
    # Import Data
    falkor = pd.read_csv("falkor_clean.csv")
    R = np.loadtxt(fpath.format("R.txt"), delimiter=",")

    # Define Inputs
    # stations = np.arange(1,15,1, dtype=float) # Stations selected for analysis
    stations = np.array((1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14))
    divider = 2  # Number of sublayers in each layer
    K = 10000  # Number of Iterations for Monte Carlo Error Propagation

    ### Data Preparation ###
    # Select Stations
    idx_station = np.where(np.isin(falkor["Station"], stations))

    # Save Data Vectors
    rho = np.array(falkor["rho"])[idx_station]  # kg/m3
    sigma0 = np.array(falkor["sigma0"])[idx_station]  # kg/m3
    DIC = np.array(falkor["DIC"])[idx_station]  # umol/kg
    DIP = np.array(falkor["DIP"])[idx_station]  # umol/kg
    NO2 = np.array(falkor["NO2"])[idx_station]  # umol/kg
    NO3 = np.array(falkor["NO3"])[idx_station]  # umol/kg
    NH4 = np.array(falkor["NH4"])[idx_station]  # umol/kg
    Nstar = np.array(falkor["Nstar"])[idx_station]  # umol/kg
    TA = np.array(falkor["TA"])[idx_station]  # umol/kg
    pH = np.array(falkor["pH"])[idx_station]
    O2 = np.array(falkor["O2"])[idx_station]  # umol/kg

    # Define Layers
    layers = np.array(
        [25, 25.5, 25.8, 26.21, 26.31, 26.44, 26.67, 27, 27.3], dtype=float
    )  # in sigma0
    sl = np.zeros((len(layers) - 1, divider + 1))
    for i in np.arange(0, len(layers) - 1):
        sl[i,] = np.linspace(layers[i], layers[i + 1], divider + 1)
    sublayers = np.unique(sl)

    ### Robust Linear Regression ###
    # Create Results Arrays
    fitting_NO3 = np.zeros((len(sublayers) - 1, 4))
    fitting_NO2 = np.zeros((len(sublayers) - 1, 4))
    fitting_NH4 = np.zeros((len(sublayers) - 1, 4))
    fitting_Nstar = np.zeros((len(sublayers) - 1, 4))
    fitting_TA = np.zeros((len(sublayers) - 1, 4))
    fitting_DIC = np.zeros((len(sublayers) - 1, 4))

    # Run Robost Linear Regression
    for i in np.arange(0, len(sublayers) - 1):
        lower_boundary = sublayers[i]
        upper_boundary = sublayers[i + 1]
        idx_layer = np.where((sigma0 > lower_boundary) & (sigma0 < upper_boundary))

        # Robust Regression
        xx = sm.add_constant(DIC[idx_layer])
        rf_NO3 = sm.RLM(NO3[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()
        rf_NO2 = sm.RLM(NO2[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()
        rf_NH4 = sm.RLM(NH4[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()
        rf_Nstar = sm.RLM(Nstar[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()
        rf_TA = sm.RLM(TA[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()
        rf_DIC = sm.RLM(DIC[idx_layer], xx, M=sm.robust.norms.HuberT()).fit()

        # Save Results
        fitting_NO3[i, :] = np.array(
            [rf_NO3.params[1], rf_NO3.bse[1], rf_NO3.params[0], rf_NO3.bse[0]]
        )
        fitting_NO2[i, :] = np.array(
            [rf_NO2.params[1], rf_NO2.bse[1], rf_NO2.params[0], rf_NO2.bse[0]]
        )
        fitting_NH4[i, :] = np.array(
            [rf_NH4.params[1], rf_NH4.bse[1], rf_NH4.params[0], rf_NH4.bse[0]]
        )
        fitting_Nstar[i, :] = np.array(
            [rf_Nstar.params[1], rf_Nstar.bse[1], rf_Nstar.params[0], rf_Nstar.bse[0]]
        )
        fitting_TA[i, :] = np.array(
            [rf_TA.params[1], rf_TA.bse[1], rf_TA.params[0], rf_TA.bse[0]]
        )
        fitting_DIC[i, :] = np.array(
            [rf_DIC.params[1], rf_DIC.bse[1], rf_DIC.params[0], rf_DIC.bse[0]]
        )

    # Save Slopes
    slopes_mean = np.array(
        [
            fitting_NO3[:, 0].T,
            fitting_NO2[:, 0].T,
            fitting_NH4[:, 0].T,
            fitting_Nstar[:, 0].T,
            fitting_TA[:, 0].T,
            fitting_DIC[:, 0].T,
        ]
    ).T
    slopes_se = np.array(
        [
            fitting_NO3[:, 1],
            fitting_NO2[:, 1],
            fitting_NH4[:, 1],
            fitting_Nstar[:, 1],
            fitting_TA[:, 1],
            fitting_DIC[:, 1],
        ]
    ).T

    ### Reaction Coefficient, Relative Contributions, and Residuals Calculations with Monte Carlo Error Propagation ###
    # Create Results Arrays
    slope_iter = np.zeros((6))
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
    slopes_obs = np.zeros((len(sublayers) - 1, 6))
    slopes_est = np.zeros((len(sublayers) - 1, 6))
    residuals = np.zeros((len(sublayers) - 1, 6))
    residuals_perc = np.zeros((len(sublayers) - 1, 6))
    nstar_slope = np.zeros((K))

    # Calculate Coeffs, Relative Contributions, and Residuals
    # do Monte Carlo simulation for each layer
    for i in np.arange(0, len(sublayers) - 1):
        # Run Monte Carlo
        for k in range(K):  # K is the number of iterations
            for j in range(6):
                # pick a value for measured Δtracer:DIC based on normal distribution
                slope_iter[j] = np.random.normal(
                    loc=slopes_mean[i, j], scale=slopes_se[i, j], size=1
                )
            # Calculate Coeffs
            c_temp, rnorm = sc.nnls(R, slope_iter)  # solve for Χ matrix
            coeff_iter[k, :] = c_temp / np.sum(c_temp)  # make sure X matrix sums to 1

            nstar_slope[k] = np.random.normal(
                loc=slopes_mean[i, 4], scale=slopes_se[i, 4], size=1
            )
            # Calculate Relative Importances
            anmx_iter[k] = (
                (R[3, 2] * c_temp[2])
                / (R[3, 2] * c_temp[2] + R[3, 1] * c_temp[1])
                * 100
            )
            denit_iter[k] = (
                (R[3, 1] * c_temp[1])
                / (R[3, 2] * c_temp[2] + R[3, 1] * c_temp[1])
                * 100
            )
            ox_iter[k] = (
                (R[1, 3] * c_temp[3] - R[0, 2] * c_temp[2])
                / (
                    R[1, 2] * c_temp[2]
                    + R[1, 1] * c_temp[1]
                    + R[1, 3] * c_temp[3]
                    - R[0, 2] * c_temp[2]
                )
                * 100
            )
            red_iter[k] = (
                (R[1, 2] * c_temp[2] + R[1, 1] * c_temp[1])
                / (
                    R[1, 2] * c_temp[2]
                    + R[1, 1] * c_temp[1]
                    + R[1, 3] * c_temp[3]
                    - R[0, 2] * c_temp[2]
                )
                * 100
            )
            nitox_iter[k] = (
                abs(R[1, 3] * c_temp[3])
                / abs(R[1, 3] * c_temp[3] - R[1, 0] * c_temp[0])
                * 100
            )
            dnrn_iter[k] = (
                abs(R[1, 0] * c_temp[0])
                / abs(R[1, 3] * c_temp[3] - R[1, 0] * c_temp[0])
                * 100
            )
            caco3_iter[k] = (
                abs(c_temp[4])
                / (abs(c_temp[0] + c_temp[1] - c_temp[2] - c_temp[3]) + c_temp[4])
                * 100
            )
            otherDIC_iter[k] = (
                abs(c_temp[0] + c_temp[1] - c_temp[2] - c_temp[3])
                / (abs(c_temp[0] + c_temp[1] - c_temp[2] - c_temp[3]) + c_temp[4])
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
                abs(R[1, 2] * c_temp[2])
                / abs(
                    -R[1, 2] * c_temp[2]
                    - R[1, 1] * c_temp[1]
                    - R[1, 3] * c_temp[3]
                    + R[1, 0] * c_temp[0]
                )
                * 100
            )
            relnit_denit_iter[k] = (
                abs(R[1, 1] * c_temp[1])
                / abs(
                    -R[1, 2] * c_temp[2]
                    - R[1, 1] * c_temp[1]
                    - R[1, 3] * c_temp[3]
                    + R[1, 0] * c_temp[0]
                )
                * 100
            )
            relnit_nitox_iter[k] = (
                abs(R[1, 3] * c_temp[3])
                / abs(
                    -R[1, 2] * c_temp[2]
                    - R[1, 1] * c_temp[1]
                    - R[1, 3] * c_temp[3]
                    + R[1, 0] * c_temp[0]
                )
                * 100
            )
            relnit_dnrn_iter[k] = (
                abs(R[1, 0] * c_temp[0])
                / abs(
                    -R[1, 2] * c_temp[2]
                    - R[1, 1] * c_temp[1]
                    - R[1, 3] * c_temp[3]
                    + R[1, 0] * c_temp[0]
                )
                * 100
            )
            relnit_iter[k, :] = [
                relnit_anmx_iter[k],
                relnit_denit_iter[k],
                relnit_nitox_iter[k],
                relnit_dnrn_iter[k],
            ]

        # Save the Results
        # take the mean of the 1000 Monte Carlo iterations
        coeff_mean[i, :] = np.mean(coeff_iter, axis=0)  
        coeff_se[i, :] = np.std(coeff_iter, axis=0)  # /math.sqrt(K)
        relimp_mean[i, :] = np.nanmean(relimp_iter, axis=0)
        relimp_se[i, :] = np.nanstd(relimp_iter, axis=0)  # /math.sqrt(K)
        relnit_mean[i, :] = np.nanmean(relnit_iter, axis=0)
        relnit_se[i, :] = np.nanstd(relnit_iter, axis=0)  # /math.sqrt(K)

        # Calculate Residuals
        slopes_obs[i, :] = slopes_mean[i, :].T
        slopes_est[i, :] = np.dot(R, coeff_mean[i, :])
        residuals[i, :] = slopes_obs[i, :] - slopes_est[i, :]
        residuals_perc[i, :] = abs(residuals[i, :]) / abs(slopes_obs[i, :]) * 100

    ### Save the Outputs as a CSV file ###
    np.savetxt(
        fpath.format("slopes_mean.csv"),
        slopes_mean,
        delimiter=",",
        header="NO3, NO2, NH4, Nstar, TA, DIC",
    )
    np.savetxt(
        fpath.format("slopes_se.csv"),
        slopes_se,
        delimiter=",",
        header="NO3, NO2, NH4, Nstar, TA, DIC",
    )
    np.savetxt(
        fpath.format("coeffs_mean.csv"),
        coeff_mean,
        delimiter=",",
        header="DNRN, Anammox, Denitrificatio, Nitrite Oxidation, CaCO3 Dissolution ",
    )
    np.savetxt(
        fpath.format("coeffs_se.csv"),
        coeff_se,
        delimiter=",",
        header="DNRN, Anammox, Denitrificatio, Nitrite Oxidation, CaCO3 Dissolution ",
    )
    np.savetxt(
        fpath.format("relative_importances_mean.csv"),
        relimp_mean,
        delimiter=",",
        header="anmx, denit, ox, red, nitox, dnrn, caco3 diss, other DIC",
    )
    np.savetxt(
        fpath.format("relative_importances_nitrite_mean.csv"),
        relnit_mean,
        delimiter=",",
        header="anmx, denit, nitox, dnrn",
    )
    np.savetxt(
        fpath.format("relative_importances_nitrite_se.csv"),
        relnit_se,
        delimiter=",",
        header="anmx, denit, nitox, dnrn",
    )
    np.savetxt(
        fpath.format("relative_importances_se.csv"),
        relimp_se,
        delimiter=",",
        header="anmx, denit, ox, red, nitox, dnrn, caco3 diss, other DIC",
    )
    np.savetxt(
        fpath.format("residuals.csv"),
        residuals,
        delimiter=",",
        header="NO3, NO2, NH4, Nstar, TA, DIC",
    )
    np.savetxt(
        fpath.format("percentage_residuals.csv"),
        residuals_perc,
        delimiter=",",
        header="NO3, NO2, NH4, Nstar, TA, DIC",
    )
    np.savetxt(
        fpath.format("reaction_matrix.csv"),
        R,
        delimiter=",",
        header="Columns: Reactions, Rows: Tracers",
    )
    np.savetxt(fpath.format("layers.csv"), sublayers, delimiter=",", header="Layers")
    print("done")


if __name__ == "__main__":
    # loop through all Cox scenarios and recalculate outputs
    ignore = [".DS_Store", "redfield", "experimental2"]
    output_dirs = [
        ele for ele in os.listdir("output/OM_variations/") if ele not in ignore
    ]
    for d in output_dirs:
        print(d)
        fpath = "output/OM_variations/" + d + "/{}"
        calc_falkor_outputs(fpath)
