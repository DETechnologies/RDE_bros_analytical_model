%% these are plotsin the paper!
close all force
clear 
clc
%  in the no-delete folder we have data from 

% %% This SECTION IS FOR VARYING MDOT ONLY!!
% % THIS PLOTS THE THRUST MDOT AND ISP FIGURE ----------------------------------------------------------------
% load('noDelete_mat\Output_data_jan04_vMdot_LP.mat')
% figure("Name","Thrust Mdot Isp @P0=1000kPa, T0=290K")
% grid on
% hold on
% plot(Output(:,37),Output(:,31));
% ylabel('Thrust [N]')
% % ylim([600 2200])
% 
% yyaxis right
% plot(Output(:,37),Output(:,35));
% ylim([428 470])
% ylabel('Specific Impulse, I_{SP} [s]')
% 
% xlabel('Mass flowrate, m_{dot} [kg/s]')
% legend('m_{dot}[kg/s]','I_{SP}[s]')
% xlim([0.18 0.4])
% 
% % This plots the comparison data on top of ours.
% x=readtable('noDelete_mat\AmericanEngineData.xlsx');
% xlsx=table2array(x(:,1:3));
% 
% AFRL=xlsx(1:20,1:3);
% Purdue=xlsx(21:27,1:3);
% UCF=xlsx(28:40,1:3);
% 
% figure("Name","Thrust")
% grid on
% hold on
% plot(Output(:,37),Output(:,31),color="black");
% ylabel('Thrust [N]')
% 
% scatter(AFRL(:,1),AFRL(:,2))
% scatter(Purdue(:,1),Purdue(:,2))
% scatter(UCF(:,1),UCF(:,2))
% 
% xlabel('Mass flowrate, m_{dot} [kg/s]')
% legend('m_{dot}[kg/s]','AFRL','Purdue','UCF')
% xlim([0.25 0.40])
% ylim([0 2000])
% 
% figure("Name","ISP")
% grid on
% hold on
% plot(Output(:,37),Output(:,35),color="black");
% hold on
% ylabel('Specific Impulse, I_{SP} [s]')
% 
% scatter(AFRL(:,1),AFRL(:,3),color="r")
% scatter(Purdue(:,1),Purdue(:,3),color="y")
% scatter(UCF(:,1),UCF(:,3),color="m")
% 
% xlabel('Mass flowrate, m_{dot} [kg/s]')
% legend('I_{SP}[s]','AFRL','Purdue','UCF')
% xlim([0.25 0.40])
% ylim([0 600])
%% THIS SECTION IS FOR VARYING EQV RATIO ONLY!!!
load('Output_combined_EqvRvary_yesSean.mat')
% THIS PLOTS THE CELL WIDTH FIGURE ----------------------------------------------------------------
figure("Name","Cell Sizes - big")
plt=loglog(Output(:,3),Output(:,16:19));
hold on
grid on

% set(gca, 'XTick', [1:10])
set(gca, 'XTick', [0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10])
xlabel("Equivalence Ratio φ")
ylabel("Detonation Cell Size, λ [m]")
legend(Output_dataNames(1,16:19))
% % do not include sean CB in this plot it just makes a straight line.
% legend("Westbrook Cell Size [m]","Cell Width Estimate (A=29)","Gavrikov Cell Size [m]","Ng. et al. Cell Size [m]", "SeanCB Cell Size [m]")
% xlim([0.1 10])

% THIS PLOTS EQV_R WITH A THRESHOLD LINE SHOWING WHERE 1.6MM IS.
figure("Name","Cell Sizes - threshold")
loglog(Output(:,3),Output(:,16),color=[0, 0.4470, 0.7410]);
loglog(Output(:,3),Output(:,16:19));

hold on
grid on
yline(1.6/1000)

xlabel("Equivalence Ratio φ")
ylabel("Detonation Cell Size, λ [m]")
% legend("Westbrook Cell Size [m]","Minimum Cell Size (1.6mm)")
legend("Westbrook Cell Size [m]","Cell Width Estimate (A=29)","Gavrikov Cell Size [m]","Ng. et al. Cell Size [m]", "Minimum Cell Size (1.6mm)")

% xlim([0.4 1.18])
% ylim([1.075/1000 1.5/1000])

% THIS PLOTS GEOMETRY RELATED THINGS in the following
% plots!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

figure("Name","Geometry (all)")
plot(Output(:,3),Output(:,20:27));

hold on
grid on
yline(1.6/1000)
legend(Output_dataNames(1,20:27))
xlim([0.9 1.1])

% Plotting the OD
figure("Name","Geometry-OD")
plot(Output(:,3),Output(:,[24,20]));
hold on
grid on
xline(1.015)
xlabel("Equivalence Ratio φ")
ylabel("Minimum Outer Diameter, OD_{min} [mm]")
legend("Minimum Channel OD (method 1) [mm]","Minimum Channel OD (method 2) [mm]","Maximum Equivalence Ratio φ (1.015)")
xlim([0.65 1.05])

% PLOTTING MINIMUM DELTA
figure("Name","Geometry-Delta")
plot(Output(:,3),Output(:,[25,21]));
hold on
grid on
xline(1.015)
xlabel("Equivalence Ratio φ")
ylabel("Minimum Chamber Width, δ [mm]")
legend("Minimum Channel Width (method 1) [mm]","Minimum Channel Width (method 2) [mm]","Maximum Equivalence Ratio φ (1.015)")
xlim([0.65 1.05])

% PLOTTING MINIMUM length
figure("Name","Geometry-length")
plot(Output(:,3),Output(:,[26,22]));
hold on
grid on
xline(1.015)
xlabel("Equivalence Ratio φ")
ylabel("Minimum Length L_{min} [mm]")
legend("Minimum Channel Length (method 1) [mm]","Minimum Channel Length (method 2) [mm]","Maximum Equivalence Ratio φ (1.015)")
xlim([0.65 1.05])

% PLOTTING MINIMUM FILL HEIGHT
figure("Name","Geometry-fill height")
plot(Output(:,3),Output(:,23));
hold on
grid on
xline(1.015)
xlabel("Equivalence Ratio φ")
ylabel("Minimum Fill Height, h^* [mm]")
% legend(Output_dataNames(1,[20,24]))
legend("Minimum Fill Height (method 1) [mm]","Maximum Equivalence Ratio φ (1.015)")
xlim([0.65 1.05])
%% vary pressure 
load('OutputVPressure_v3_combined.mat')

% this is wavespeed and mdot plot in the paper ----------------------------------------------------------------
figure("Name","Wavespeed and M_dot, T0=290K")
grid on
hold on
plot(Output(:,30),Output(:,6),Marker="x");
ylabel('CJ Wave Speed [m/s]')
% ylim([ 1800])
xlabel('Mass flowrate, m_{dot} [kg/s]')
% xlim([0.22 0.4])

%THIS PLOTS THE CELL WIDTH FIGURE ----------------------------------------------------------------
figure("Name","Cell Sizes varying pressure")
plt=plot(Output(:,1),Output(:,16:19));
grid on

xlim([0.75*101.325e+3 2*101.325e+3])
xlabel("Input Pressure P_0 [Pa]")
ylabel("Detonation Cell Size, λ [m]")
legend(Output_dataNames(1,16:19))

% this plots the thrust and pressure
figure("Name","pressure vs thrust")
plot(Output(:,1),Output(:,31));
hold on
plot([24e+3;880e+3],[1350;1350],Color="black"); % plot horizontal line at 1350N
xlim([24e+3 880e+3])
xlabel("Stagnation Pressure P_0 [Pa]")
ylabel("Thrust [N]")
% legend("Thrust (m_{dot} into engine)","Thrust Objective (1350N)")

plot([101.325e+3;101.325e+3],[1000;2500],Color="black",LineStyle="--"); %this plots a vertical line showing where atmostpheric pressure is.
yyaxis right
ylabel("Detonation Cell Size \lambda [mm]")

loglog(Output(:,1),Output(:,37),Color="black",LineStyle=":",LineWidth=2);
legend("Thrust (m_{dot} into engine)","Thrust Objective (1350N)","Atmospheric Pressure (101.325kPa)","SeanCB Cell Size")


%% this section plots things with varying input temperatures
load('OutputVtemp_try2_YesSean.mat')
figure("Name","Cell Sizes varying temp")
plt=loglog(Output(:,2),Output(:,16:19));
grid on

xlabel("Input Temperature T_0 [°K]")
ylabel("Detonation Cell Size, λ [m]")
legend(Output_dataNames(1,16:19))
xlim([273 320])
