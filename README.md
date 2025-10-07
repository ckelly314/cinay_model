1. Run ```formatdata.py```
	Read in updated results ("5906484qcno2_updated.txt"), reformat, and calculate inputs for CANYON-B algorithm.

2. Copy input_for_CANYONB.csv over to CANYON-B Argo folder.

3. Run ```run_CANYON_Argo.m``` (in CANYON-B Argo folder)
	Run CANYON-B neural network on Argo data to predict phosphate.

4. Copy output ("CANYONB_output.csv") back over to Cinay_model_updated folder.

5. Run ```clean_argo.py```
    Renames columns, selects data flagged as good, and saves the analyzed variables and data as a .csv file.

6. Run ```calc_oxycline_features_argo.ipynb```
    Estimates organic matter C:N:P ratio and AOU in the oxycline using a robust least squares model in the oxycline (depth < 200 and 2 < [O2] < 190).

7. Run ```calc_reaction_stoi_and_R.py```
    Calculates and saves reaction stoichiometric coefficients and R matrices for different organic matter compositions, including the "experimental" C:N:P determined in calc_oxycline_features_argo.ipynb. The R matrix is the expected change in (NO3, NO2, NH4, N*, TA, and DIC) normalized to the change in DIC for (NO3- -> NO2-, NO2- -> N2O, Anammox, NO2- -> NO3-, and CaCO3 dissolution) (Table 1 of Cinay et al., 2022).

8. Run ```calc_argo_outputs.py```
	Calculates the relative reaction rates and relative contributions to CaCO3 dissolution, based on the R matrix output by calc_reaction_stoi_and_R.py. Saves relative contributions as a text file for further plotting.

9. Run ```plot_argo_results_Cox.py```
	Plots (i) the relative reaction rates (coefficients, i.e. Chi matrix in paper) as a heat map, (ii) residuals compared to the layer by layer tracer slope data, and (iii) relative contributions for a given output from calc_falkor_outputs.py.
