%% these are plotsin the paper!
close all force
clear 
clc
%  in the no-delete folder we have data from 

%% This is the thrust and isp versus mdot plot in the paper
% THIS PLOTS THE THRUST MDOT AND ISP FIGURE ----------------------------------------------------------------
load('noDelete_mat\Output_data_dec17_veqvR.mat')
figure("Name","Thrust Mdot Isp @P0=1000kPa, T0=290K")
grid on
hold on
plot(Output(:,30),Output(:,31),Marker="x");
ylabel('Thrust [N]')
ylim([400 1800])

yyaxis right
plot(Output(:,30),Output(:,35),Marker="o");
ylim([425 490])
ylabel('Specific Impulse, I_{SP} [s]')

xlabel('Mass flowrate, m_{dot} [kg/s]')
legend('m_{dot}[kg/s]','I_{SP}[s]')
xlim([0.18 0.4])

% THIS PLOTS THE CELL WIDTH FIGURE ----------------------------------------------------------------
figure("Name","Cell Sizes")
plt=loglog(Output(:,3),Output(:,16:19),LineStyle="--",Marker="o");
grid on

xlabel("Equivalence Ratio")
ylabel("Cell Size [m]")
legend(Output_dataNames(1,16:19))

% THIS PLOTS THE CELL WIDTH FIGURE BETTER RESOLUTION  ----------------------------------------------------------------
figure("Name","Cell Sizes")
plt=loglog(Output(:,3),Output(:,16:19),LineStyle="--",Marker="o");
grid on

xlabel("Equivalence Ratio")
ylabel("Cell Size [m]")
legend(Output_dataNames(1,16:19))
xlim([0.6 1.5])

%% vary pressure 
load('noDelete_mat\Output_data_dec14_vpressure.mat')

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
figure("Name","Cell Sizes")
plt=loglog(Output(:,1),Output(:,16:19),LineStyle="--",Marker="o");
grid on

xlabel(Output_dataNames(1,1))
ylabel("Cell Size [m]")
legend(Output_dataNames(1,16:19))




%% this section plots things with varying input temperatures
load('noDelete_mat\Output_data_dec14_vpressure.mat')
