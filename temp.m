%% temp test delete after
clear 
clc 
close all force

% vary pressure 
load('noDelete_mat\Output_data_dec17_vpressure.mat')


%THIS PLOTS THE CELL WIDTH FIGURE ----------------------------------------------------------------
figure("Name","V PRESSURE LOGLOG")
plt=loglog(Output(:,1),Output(:,16:19),LineStyle="--",Marker="o");
grid on

xlabel("Input Pressure P_0 [Pa]")
ylabel("Detonation Cell Size, λ [m]")
legend(Output_dataNames(1,16:19))


figure("Name","V PRESSURE REGULA")
plt=plot(Output(:,1),Output(:,16:19),LineStyle="--",Marker="o");
grid on

xlabel("Input Pressure P_0 [Pa]")
ylabel("Detonation Cell Size, λ [m]")
legend(Output_dataNames(1,16:19))