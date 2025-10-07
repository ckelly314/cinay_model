% add path to CO2SYS
addpath('/Users/colette/Drive/Shared Drives/GO-BGC-Workshop-2023/GO-BGC-2023 Workshop/Project Teams/Group3_Bourbonnais_Bif_Altabet/2024_ODZ_Paper/Model/CO2SYS-MATLAB-v1/src')

data = readtable('input_for_CANYONB.csv');

out = CANYONB(data.gtime,data.lat,data.lon,data.pres,data.temp,data.psal,data.doxy,{'NO3','PO4','SiOH4'});

data.PO4_CANY = out.PO4;
data.NO3_CANY = out.NO3;
data.SiOH4_CANY = out.SiOH4;

writetable(data, 'CANYONB_output.csv');
