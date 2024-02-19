%% Compare detonation cell size predictions (westbrook) from detonation database (caltech)
% to the detonation cell size predictions from this model.

close all force
clear
clc

%% load data
load("Output_data_feb19_varyPressure.mat");
detonationDatabase=readtable("CellSizes_combined_vs_initPressure.xlsx");


%% plot the things on a loglog plot

figure("Name","Detonation Cell Sizes versus Init Pressure")
scatter(table2array(detonationDatabase(:,1)),table2array(detonationDatabase(:,2)));
hold on
scatter((Output(:,1)/1000),(Output(:,16)*1000))
grid on
set(gca,'xscale','log','yscale','log')

xlabel('initialPressure [kpa]')
ylabel('Cell size [mm]')

legend("DetonationDatabase","Calculator - Westbrook")