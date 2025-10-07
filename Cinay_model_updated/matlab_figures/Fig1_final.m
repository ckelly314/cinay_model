%% Fig.1 Map for nitrite, highlighting the cruise and float tracking - MBB 03/14/2025
% This code uses cruise tracking and float tracking to map their locations.
% We then calculate poc and chla stocks and fit fig.1b and fig.1c afterwards.


close all
clear all
clc


%% Define variables for the map visualization (Fig.1a)
% Open excel spreadsheet (Altabet SR2114 cruise) and textfile (float 5906484)

% Cruise:

cruisetable = 'Altabet SR2114 NO3 and NO2 data.xlsx';
cruisetable = readtable(cruisetable);
% Extract columns from cruise
Cstation = cruisetable.Var2;
Clatitude = cruisetable.Var4;
Clongitude = cruisetable.Var5;

% Find unique stations
[uniqueStations, ia] = unique(Cstation);
% Extract corresponding latitude and longitude
uniqueLatitudes = Clatitude(ia);
uniqueLongitudes = Clongitude(ia);
% Create a table for cruise coordinates
Cposition = table(uniqueStations, uniqueLatitudes, uniqueLongitudes, ...
                  'VariableNames', {'Station', 'Latitude', 'Longitude'});
Cposition.Longitude = mod(360 + Cposition.Longitude, 360); % conversion from degrees west to east


% Float:

float = '5906484qcno2.txt';
float = readtable(float, 'Delimiter', '\t'); 
% Extract columns from float
Station = float.Station;
Lat = float.Lat__N_;
Lon = float.Lon__E_;
% Find unique stations
[uniqueStations, ia] = unique(Station);
% Extract corresponding latitude and longitude
uniqueLatitudes = Lat(ia);
uniqueLongitudes = Lon(ia);
% Create a table with float coordinates
Fposition = table(uniqueStations, uniqueLatitudes, uniqueLongitudes, ...
                  'VariableNames', {'Station', 'Latitude', 'Longitude'});

% Determine the range of Latitude and Longitude for map
lat_range = [5 20];
lon_range = [min(Cposition.Longitude)+5 max(Cposition.Longitude)+2];

% Create subplots with specified layout
figure;

% Upper subplot for the map (larger)
subplot(3, 1, 1);
m_proj('Mercator', 'lat', lat_range, 'lon', lon_range);
m_grid('box', 'fancy', 'tickdir', 'out');
hold on;
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'none');

% Remove the 5 highest latitude values from Cposition
[~, idx] = sort(Cposition.Latitude, 'descend');
Cposition(idx(1:5), :) = [];

% Plot
h3 = m_scatter(Fposition.Longitude, Fposition.Latitude, 20, 'r', 'filled');
h1 = m_scatter(Cposition.Longitude, Cposition.Latitude, 20, 'k','LineWidth', 2);
legend([h1, h3], {'SR2114 Cruise', 'Float 5906484'}, 'Location', 'northeast');
text(0.02, 0.3, 'a', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

%% POC and Chla stocks over time and seasons (Fig.1b and Fig.1c)

data = readtable('5906484qcno2.txt', 'Delimiter', '\t');
data.mon_day_yr = datetime(data.mon_day_yr, 'InputFormat', 'MM/dd/yyyy');
unique_dates = unique(data.mon_day_yr);

% Initialize arrays to store integrated Chl_a and POC values and corresponding dates
integrated_Chl_a = [];
integrated_POC = [];
date_list = [];

for i = 1:length(unique_dates)
    current_date = unique_dates(i);
    date_data = data(data.mon_day_yr == current_date, :);
    
    % Check if date_data is empty
    if isempty(date_data)
        fprintf('No data found for date: %s\n', char(current_date));
        continue;
    end
    
    % Filter the data for depths less than or equal to 150 meters
    upper_150m_data = date_data(date_data.Depth_m_ <= 150, :);
   
    if isempty(upper_150m_data)
        fprintf('No data found for depths <= 150 meters on date: %s\n', char(current_date));
        continue;
    end
    upper_150m_data = sortrows(upper_150m_data, 'Depth_m_');
    
    % Integrate the Chl_a values for the upper 150 meters
    integrated_value_Chl_a = trapz(upper_150m_data.Depth_m_, upper_150m_data.Chl_a_mg_m_3_);
    integrated_Chl_a = [integrated_Chl_a; integrated_value_Chl_a];
    
    % Integrate the POC values for the upper 150 meters
    integrated_value_POC = trapz(upper_150m_data.Depth_m_, upper_150m_data.POC_mmol_m_3_);
    integrated_POC = [integrated_POC; integrated_value_POC];
    
    date_list = [date_list; current_date];
end

% Determine the color for each date based on the season
colors = zeros(length(date_list), 3);
for i = 1:length(date_list)
    month_value = month(date_list(i));
    if ismember(month_value, [3, 4, 5])
        colors(i, :) = [0, 1, 0]; % Green for Spring
    elseif ismember(month_value, [6, 7, 8])
        colors(i, :) = [1, 0.5, 0]; % Orange for Summer
    elseif ismember(month_value, [9, 10, 11])
        colors(i, :) = [0.5, 0, 0.5]; % Purple for Fall
    else
        colors(i, :) = [0, 0, 1]; % Blue for Winter
    end
end

% Calculate moving average
window_size = 6;
movmean_values_Chl_a = movmean(integrated_Chl_a, window_size);
movmean_values_POC = movmean(integrated_POC, window_size);

% Middle subplot for Chl_a 
subplot(3, 1, 2);
yyaxis left
b1 = bar(date_list, integrated_Chl_a, 'FaceColor', 'flat');
b1.CData = colors;
ylabel({'Chla in the upper'; '150 m (mg m^-^2)'});

yyaxis right
plot(date_list, movmean_values_Chl_a, '-k','LineWidth', 2);
ylabel({'Moving Average'; '(mg m^-^2)'});
ax1 = gca;
ax1.YColor = 'k';

text(0.02, 0.95, 'b', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

xtickformat('MM/dd/yy');
xtickangle(45);

% Add legend for seasons
hold on;
spring_patch = patch([NaN NaN], [NaN NaN], [0, 1, 0], 'DisplayName', 'Spring');
summer_patch = patch([NaN NaN], [NaN NaN], [1, 0.5, 0], 'DisplayName', 'Summer');
fall_patch = patch([NaN NaN], [NaN NaN], [0.5, 0, 0.5], 'DisplayName', 'Fall');
winter_patch = patch([NaN NaN], [NaN NaN], [0, 0, 1], 'DisplayName', 'Winter');
legend([spring_patch, summer_patch, fall_patch, winter_patch], 'Location', 'northeast');
hold off;

% Lower subplot for POC 
subplot(3, 1, 3);
yyaxis left
b2 = bar(date_list, integrated_POC, 'FaceColor', 'flat');
b2.CData = colors;
ylim([0 500])
ylabel({'POC in the upper'; '150 m (mmol m^-^2)'});

yyaxis right
plot(date_list, movmean_values_POC, '-k', 'LineWidth', 2);
ylabel({'Moving Average'; '(mmol m^-^2)'});
ax2 = gca;
ax2.YColor = 'k'; 

xtickformat('MM/dd/yy');
xtickangle(45);

text(0.02, 0.95, 'c', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

saveas(gcf,'Fig_1','tiff'); saveas(gcf,'Fig_1','svg'); saveas(gcf,'Fig_1.fig');
