%% Compare detonation cell size predictions (westbrook) from detonation database (caltech)
% to the detonation cell size predictions from this model.

close all force
clear
clc

%% load data
% load("Output_data_feb23_combined.mat");
load('..\noDelete_mat\OutputVpressure_try2_yesSean2.mat')
detonationDatabase=readtable("CellSizes_combined_vs_initPressure.xlsx");


%% plot the things on a loglog plot

figure("Name","Detonation Cell Sizes versus Init Pressure")
scatter(table2array(detonationDatabase(:,1)),table2array(detonationDatabase(:,2)),'x','black');
hold on
loglog((Output(:,1)/1000),(Output(:,16)*1000),LineWidth=2)%,Color=[.5 .5 .5])
% loglog((Output(:,1)/1000),(Output(:,17:19)*1000),LineWidth=2)
loglog((Output(:,1)/1000),(Output(:,18:19)*1000),LineWidth=2)
loglog((Output(:,1)/1000),(Output(:,37)),LineWidth=2)
grid on
set(gca,'xscale','log','yscale','log')

xlabel('Stagnation Pressure, P_0 [kpa]')
ylabel('Cell Size \lambda [mm]')
xlim([4 1215])

% legend("CalTech Detonation Database","Westbrook","A=29","Garikov","Ng et al.","Sean CB")
legend("CalTech Detonation Database","Westbrook","Garikov","Ng et al.","Sean CB")










