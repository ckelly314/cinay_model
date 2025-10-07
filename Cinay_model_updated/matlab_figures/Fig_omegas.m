%% Calcite and Aragonite Omegas - MBB 041625

% Estimating omega for the time-series using CO2SYS
% Sharp, J. D., Pierrot, D., Humphreys, M. P., Epitalon, J.-M., Orr, J. C., Lewis, E. R., 
% & Wallace, D. W. R. (2023). CO2SYSv3 for MATLAB (v3.2.1). Zenodo. https://doi.org/10.5281/zenodo.7552554

close all
clear all
clc

data = readtable('5906484qcno2.txt', 'Delimiter', '\t');

% Filtering bad pH
station_56_idx = data.Station == 56;
station_57_idx = data.Station == 57;
pH_station_56 = data.pHinsitu_Total_(station_56_idx);
data.pHinsitu_Total_(station_57_idx) = pH_station_56;

%% Parameters for CO2SYS
PAR1 = data.TALK_LIAR__mol_kg_;
PAR2 = data.pHinsitu_Total_;
SAL = data.Salinity_pss_;
TEMPIN = data.Temperature__C_;
PRESIN = 0; % Surface pressure 
PAR1TYPE = 1; % ALKALINITY type
PAR2TYPE = 3; % pH type
pHSCALEIN = 1;
K1K2CONSTANTS = 17;
KSO4CONSTANT = 1;
KFCONSTANT = 2;
BORON = 2;

% Calculate Omegas: Columns 17 (Calcite) and 18 (Aragonite)
[DATA_CO2SYS,HEADERS,NICEHEADERS] = CO2SYS_adjusted_to_v2_0_5(PAR1, PAR2, PAR1TYPE,...
    PAR2TYPE, SAL, TEMPIN, nan, PRESIN, nan, 1, 1, 0,0, pHSCALEIN, K1K2CONSTANTS,...
    KSO4CONSTANT, KFCONSTANT, BORON);

%% Plotting omegas vs depth and time-series

Depth_m = data.Depth_m_;
sigma_prof = data.Sigma_theta_kg_m_3_;
datet = data.mon_day_yr;
omegaCa = DATA_CO2SYS(:, 17);
omegaAr = DATA_CO2SYS(:, 18);

% Substitute bad values with NaN
omegaCa(omegaCa == -999) = NaN;
omegaAr(omegaAr == -999) = NaN;

%% Find shallowest depth where omegas cross or approach 1
unique_stations = unique(data.Station);

depth_omegaCa_1 = NaN(size(unique_stations));
depth_omegaAr_1 = NaN(size(unique_stations));
sigma_prof_omegaCa_1 = NaN(size(unique_stations));
sigma_prof_omegaAr_1 = NaN(size(unique_stations));

for i = 1:length(unique_stations)
    station_idx = data.Station == unique_stations(i);

    % Grab station-specific data
    depth_station = Depth_m(station_idx);
    omegaCa_station = omegaCa(station_idx);
    omegaAr_station = omegaAr(station_idx);
    sigma_station = sigma_prof(station_idx);

    % Sort by increasing depth
    [depth_sorted, sort_idx] = sort(depth_station);
    omegaCa_sorted = omegaCa_station(sort_idx);
    omegaAr_sorted = omegaAr_station(sort_idx);
    sigma_sorted = sigma_station(sort_idx);

    % Find first depth where omega <= 1 for Ca
    idx_Ca = find(omegaCa_sorted <= 1, 1, 'first');
    if isempty(idx_Ca)
        [~, idx_Ca] = min(abs(omegaCa_sorted - 1));
    end
    depth_omegaCa_1(i) = depth_sorted(idx_Ca);
    sigma_prof_omegaCa_1(i) = sigma_sorted(idx_Ca);

    % Find first depth where omega <= 1 for Ar
    idx_Ar = find(omegaAr_sorted <= 1, 1, 'first');
    if isempty(idx_Ar)
        [~, idx_Ar] = min(abs(omegaAr_sorted - 1));
    end
    depth_omegaAr_1(i) = depth_sorted(idx_Ar);
    sigma_prof_omegaAr_1(i) = sigma_sorted(idx_Ar);
end

unique_dates = unique(datet);

%% Plotting

% Define productive period and nitrite regime switch
start_date = datetime('07/30/2022','InputFormat','MM/dd/yyyy');
end_date = datetime('10/11/2022','InputFormat','MM/dd/yyyy');
switch_date = datetime('11/11/2022','InputFormat','MM/dd/yyyy');

figure;
tiledlayout(2, 2);

% Left side plot
nexttile([2, 1]);
hold on;
scatter(omegaCa, Depth_m, 'b^', 'SizeData', 20); 
scatter(omegaAr, Depth_m, 'r^', 'SizeData', 20); 
set(gca, 'YDir', 'reverse');
set(gca, 'XAxisLocation', 'top');
xline(1, 'k-', 'LineWidth', 2);
xlim([0.5 1.5]);
xticks(0.5:0.2:1.5);
ylim([0 800]);
ylabel('Depth (m)');
title('\color{black}\Omega \color{blue}Calcite \color{black}and \color{red}Aragonite', 'FontSize', 14);
set(gca, 'LineWidth', 1); 
hold off;

% Upper plot on the right side
nexttile;
hold on;

% Shaded productive period
ylims = [min(depth_omegaCa_1(:)) - 20, max(depth_omegaCa_1(:)) + 20]; % a little padding
fill([start_date end_date end_date start_date], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.8 0.8 0.8], ...
    'EdgeColor', 'k', 'LineStyle', '--', 'FaceAlpha', 0.3);

% Nitrite regime division line
xline(switch_date, 'Color', [0.65 0.16 0.16], 'LineStyle', '--', 'LineWidth', 2);


yyaxis left;
plot(unique_dates, depth_omegaCa_1, 'b^-', 'LineWidth', 1.5); % Blue line
set(gca, 'YDir', 'reverse');
ylabel('Depth (m) where \Omega_C_a \approx 1', 'Color', 'b');
ax = gca;
ax.YColor = 'b';

yyaxis right;
plot(unique_dates, sigma_prof_omegaCa_1, 'k-', 'LineWidth', 1.5); 
set(gca, 'YDir', 'reverse');
ylabel('\sigma_{\theta} (kg m^{-3})', 'Color', 'k');
ax = gca;
ax.YColor = 'k';

set(gca, 'LineWidth', 1);
xtickformat('MM/yy');
xtickangle(45);
hold off;

% Bottom right: depth of omega=1 Ar vs time
nexttile;
hold on;
% Shaded productive period
ylims = [min(depth_omegaAr_1(:)) - 20, max(depth_omegaAr_1(:)) + 20]; % a little padding
fill([start_date end_date end_date start_date], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.8 0.8 0.8], ...
    'EdgeColor', 'k', 'LineStyle', '--', 'FaceAlpha', 0.3);

% Nitrite regime division line
xline(switch_date, 'Color', [0.65 0.16 0.16], 'LineStyle', '--', 'LineWidth', 2);


yyaxis left;
plot(unique_dates, depth_omegaAr_1, 'r^-', 'LineWidth', 1.5); % Red
set(gca, 'YDir', 'reverse');
ylabel('Depth (m) where \Omega_{Ar} \approx 1', 'Color', 'r');
ax = gca;
ax.YColor = 'r';

yyaxis right;
plot(unique_dates, sigma_prof_omegaAr_1, 'k-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');
ylabel('\sigma_{\theta} (kg m^{-3})', 'Color', 'k');
ax.YColor = 'k';

set(gca, 'LineWidth', 1);
xtickformat('MM/yy');
xtickangle(45);
hold off;

set(findall(gcf,'-property','FontSize'),'FontSize',14)

 saveas(gcf,'Fig_S_omegas','tiff');  saveas(gcf,'Fig_S_omegas','svg');  saveas(gcf,'Fig_S_omegas.fig');