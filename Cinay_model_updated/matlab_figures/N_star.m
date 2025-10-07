%% Time-series of N* [µmol/kg] for float WMO#5906484 - MBB 080525
% Interpolates missing N* values from bad stations (56 & 57) over time.

close all
clear all
clc

data = readtable('5906484qcno2PO4.txt', 'Delimiter', '\t');

nstar_var = 'N_star__mol_kg_';

% === Mark N* from bad stations as NaN ===
bad_nstar_idx = data.Station == 56 | data.Station == 57;
data{bad_nstar_idx, nstar_var} = NaN;

% === Create time and depth grid ===
depth_grid = 0:5:800;
time_grid = unique(data.mon_day_yr);
[T, Z] = meshgrid(datenum(time_grid), depth_grid);

% === Interpolate N* onto regular grid ===
nstar_matrix = NaN(size(T));
for i = 1:length(time_grid)
    idx = data.mon_day_yr == time_grid(i);
    z = data.Depth_m_(idx);
    nstar = data{idx, nstar_var};

    valid = ~isnan(z) & ~isnan(nstar);
    if sum(valid) > 2
        nstar_matrix(:, i) = interp1(z(valid), nstar(valid), depth_grid, 'linear', NaN);
    end
end

% === Fill missing time steps by interpolating across time ===
for k = 1:size(nstar_matrix, 1) 
    row = nstar_matrix(k, :);
    if sum(~isnan(row)) >= 3
        nstar_matrix(k, :) = interp1(find(~isnan(row)), row(~isnan(row)), ...
                                     1:length(row), 'linear', 'extrap');
    end
end

% === Plot ===
figure('Position', [100, 100, 800, 250])
pcolor(T, Z, nstar_matrix)
shading interp
set(gca, 'YDir', 'reverse')
ylabel('Depth (m)')
datetick('x', 'mm/yy', 'keeplimits')
xtickangle(45)

colormap(brewermap(16, '*YlOrRd')) 
caxis([-22 0]) 
cb = colorbar;
cb.Label.String = 'N* (µmol kg^-^1)';

hold on
[C, h] = contour(T, Z, nstar_matrix, [-20 -20], 'k--', 'LineWidth', 1.2);
hold off 

ax = gca;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
set(findall(gcf,'-property','FontSize'),'FontSize',14)

saveas(gcf, 'Fig_Nstar', 'tiff'); saveas(gcf, 'Fig_Nstar', 'svg'); saveas(gcf, 'Fig_Nstar.fig')


%% Nitrate depletion versus N*:
nitrate_var = 'Nitrate__mol_kg_';  

odz_idx = data.Depth_m_ >= 180 & data.Depth_m_ <= 300;
nstar_odz = data{odz_idx, nstar_var};

% Summary stats
fprintf('N* values between 180–300 m:\n');
fprintf('n = %d\n', numel(nstar_odz));
fprintf('Min: %.2f µmol/kg\n', min(nstar_odz, [], 'omitnan'));
fprintf('Mean: %.2f µmol/kg\n', mean(nstar_odz, 'omitnan'));
fprintf('Max: %.2f µmol/kg\n', max(nstar_odz, [], 'omitnan'));
