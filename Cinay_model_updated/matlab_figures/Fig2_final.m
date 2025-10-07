% Time-series of float WMO#5906484 - MBB 080525
% Using .txt file to plot time series with overlain potential density
% anomalies, mixed layers and ODZ.

close all
clear all
clc

data = readtable('5906484qcno2.txt', 'Delimiter', '\t');

% Filtering one bad pH. Substituting with previous profile value for
% visualization-only
% Find rows where Station is 56 and 57
station_56_idx = data.Station == 56;
station_57_idx = data.Station == 57;

pH_station_56 = data.pHinsitu_Total_(station_56_idx);% Get pHinsitu_Total_ values for Station 56
data.pHinsitu_Total_(station_57_idx) = pH_station_56;% Assign those values to Station 57

%% Estimating MLD (overlain POC)

% Extract the relevant columns
Station = data.Station;
Temperature_C = data.Temperature__C_;
Salinity_pss = data.Salinity_pss_;
Depth_m = data.Depth_m_;

% Define the temperature difference criterion - MLD
dT = 0.2;

unique_stations = unique(Station);
mld = NaN(length(unique_stations), 1);
% Calculate the MLD for each profile and potential densities
for i = 1:length(unique_stations)
    station_idx = Station == unique_stations(i);
    salt = Salinity_pss(station_idx);
    temp = Temperature_C(station_idx);
    Z = Depth_m(station_idx);
    
    % Sort the data by depth
    [Z, sortIdx] = sort(Z);
    salt = salt(sortIdx);
    temp = temp(sortIdx);
    
    % Exclude the upper 10 meters of depth
    valid_idx = Z > 10;
    salt = salt(valid_idx);
    temp = temp(valid_idx);
    Z = Z(valid_idx);
    
    mld(i) = ra_mld(salt, temp, Z, dT);
end
%% Interpolated time-series plots of float WMO#5906484 

% Variables to plot
var_names = {'Chl_a_mg_m_3_', 'POC_mmol_m_3_', 'Oxygen__mol_kg_', ...
             'pHinsitu_Total_', 'Nitrite__mol_kg_', 'Nitrate__mol_kg_'};
var_labels = {'Chl (mg m^{-3})', 'POC (mmol m^{-3})', 'O_2 (\mumol kg^{-1})', ...
              'pH', 'NO_2^- (\mumol kg^{-1})', 'NO_3 (\mumol kg^{-1})'};
ylims = [200, 800, 800, 800, 800, 800];
panel_labels = {'a', 'b', 'c', 'd', 'e', 'f'};

% Regular time and depth grid
time_grid = unique(data.mon_day_yr);
depth_grid = 0:5:800;
[T, Z] = meshgrid(datenum(time_grid), depth_grid);

% Colormaps
colormaps = {'BuGn', '*Spectral', 'PuBuGn', 'RdYlBu', '*RdYlBu', '*RdYlBu'};
% caxis ranges
caxis_ranges = {[0 0.8], [1 3], [0 3], [7.5 7.57], [0 4], [0 45]};


figure('Position', [10, 10, 800, 400]);
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
sx = gobjects(6, 1);

for v = 1:length(var_names)
    sx(v) = nexttile;
    
    var = data.(var_names{v});
    interp_matrix = NaN(size(T));
    
    % Interpolate each profile onto the depth grid
    for i = 1:length(time_grid)
        idx = data.mon_day_yr == time_grid(i);
        depth_profile = data.Depth_m_(idx);
        var_profile = var(idx);
        
        valid = ~isnan(depth_profile) & ~isnan(var_profile);
        if sum(valid) > 2
            interp_matrix(:, i) = interp1(depth_profile(valid), ...
                                           var_profile(valid), ...
                                           depth_grid, 'linear', NaN);
        end
    end
    
    % Plot
    pcolor(T, Z, interp_matrix)
    shading interp
    set(gca, 'YDir', 'reverse')
    ylim([0, ylims(v)])
    datetick('x', 'mm/yy', 'keeplimits')
    xtickangle(45);
    ylabel('Depth (m)')
    
    % Add ticks in + out, and thicker
    ax = gca;
    ax.LineWidth = 1;
    ax.TickDir = 'both';
    ax.TickLength = [0.01 0.01];

    % X-ticks more frequent
    xticks = linspace(min(datenum(time_grid)), max(datenum(time_grid)), 10);
    set(gca, 'XTick', xticks)
    set(gca, 'XTickLabel', datestr(xticks, 'mm/yy'))
    
    cb = colorbar;
    cb.Label.String = var_labels{v};
    caxis(caxis_ranges{v});
    
    % Adjust colormap
    if v == 6
        colormap(sx(v), brewermap(32, colormaps{v}));
    elseif v == 5
        original_cmap = brewermap(16, colormaps{v});
        new_cmap = original_cmap(3:end, :); 
        colormap(sx(v), new_cmap);
    else
        colormap(sx(v), brewermap(16, colormaps{v}));
    end
    
    text(max(datenum(time_grid)) - 15, ylims(v) - 30, panel_labels{v}, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontSize', 10, 'BackgroundColor', 'w', 'Margin', 1)

    % === OVERLAYS ===
    if v == 1 
        hold on
        % Sigma-theta contours
        [T_contour, Z_contour] = meshgrid(datenum(time_grid), depth_grid);
        sigma_matrix = NaN(size(T_contour));
        for i = 1:length(time_grid)
            idx = data.mon_day_yr == time_grid(i);
            depth_prof = data.Depth_m_(idx);
            sigma_prof = data.Sigma_theta_kg_m_3_(idx);
            valid = ~isnan(depth_prof) & ~isnan(sigma_prof);
            if sum(valid) > 2
                sigma_matrix(:, i) = interp1(depth_prof(valid), ...
                                              sigma_prof(valid), ...
                                              depth_grid, 'linear', NaN);
            end
        end
        % Draw sigma-theta contours
        [C1, h1] = contour(T_contour, Z_contour, sigma_matrix, [22.5 22.5], 'k');
        set(h1, 'LineWidth', 1.2, 'LineStyle', '-');
        [C3, h3] = contour(T_contour, Z_contour, sigma_matrix, [24 24], 'k');
        set(h3, 'LineWidth', 1.2, 'LineStyle', '-');

        % Plot MLD
        plot(datenum(time_grid), mld, 'r--', 'LineWidth', 1.2)

        hold off
    
    elseif v == 3 || v == 5 || v == 6 
        hold on
        % Oxygen contour overlays at O2 = 1
        [T_contour, Z_contour] = meshgrid(datenum(time_grid), depth_grid);
        oxygen_matrix = NaN(size(T_contour));
        for i = 1:length(time_grid)
            idx = data.mon_day_yr == time_grid(i);
            depth_prof = data.Depth_m_(idx);
            oxygen_prof = data.Oxygen__mol_kg_(idx);
            valid = ~isnan(depth_prof) & ~isnan(oxygen_prof);
            if sum(valid) > 2
                oxygen_matrix(:, i) = interp1(depth_prof(valid), ...
                                              oxygen_prof(valid), ...
                                              depth_grid, 'linear', NaN);
            end
        end
        [C_oxygen, h_oxygen] = contour(T_contour, Z_contour, oxygen_matrix, [1 1], 'm');
        set(h_oxygen, 'LineWidth', 1.5, 'LineStyle', '-');
        
        hold off
    end
end

 saveas(gcf,'Fig_2','tiff'); saveas(gcf,'Fig_2','svg'); saveas(gcf,'Fig_2.fig');
