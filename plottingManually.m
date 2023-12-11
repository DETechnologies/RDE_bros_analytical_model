close all force
clear 
clc

load("Output_data_varyEqvRatio.mat")

%% this plots a cute little loglog plot of cell sizes vs eqv ratio!
% figure("Name","Cell Sizes: P0=5atm, T0=313.15K")
% plt=loglog(Output(:,3),Output(:,16:19),LineStyle="--",Marker="o");
% grid on
% 
% xlabel("Equivalence Ratio")
% ylabel("Cell Size [m]")
% legend(Output_dataNames(1,16:19))

%% cj and vn pressures
% figure("Name","CJ and VN Pressures: P0=5atm, T0=313.15K")
% grid on
% hold on
% scatter(Output(:,3),Output(:,[7,9]));
% 
% xlabel(Output_dataNames(1,3))
% ylabel(Output_dataNames(1,[7,9]))
% legend(Output_dataNames(1,[7,9]))

%% CJspeed and speed of sound in that gas
figure("Name","CJ speed and Speed of Sound @P0=5atm, T0=313.15K")
grid on
hold on
scatter(Output(:,3),Output(:,[5,6]));

xlabel(Output_dataNames(1,3))
ylabel(Output_dataNames(1,[5,6]))
legend(Output_dataNames(1,[5,6]))



%% outrage
% figure("Name","something outrageous @ P0=5atm, T0=313.15K")
% grid on
% hold on
% scatter(Output(:,3),Output(:,[4:6,8,10]));
% 
% xlabel(Output_dataNames(1,3))
% % ylabel(Output_dataNames(1,3:24))
% legend(Output_dataNames(1,[4:6,8,10]))
% 
% figure("Name","something outrageous pt2 @ P0=5atm, T0=313.15K")
% grid on
% hold on
% scatter(Output(:,3),Output(:,[10:12]));
% 
% xlabel(Output_dataNames(1,3))
% % ylabel(Output_dataNames(1,3:24))
% legend(Output_dataNames(1,[10:12]))

%% dependent on shape of det cell size so it shows nothign
% figure("Name","OD: P0=5atm, T0=313.15K")
% plt=loglog(Output(:,3),Output(:,22),LineStyle="--",Marker="o");
% grid on
% 
% xlabel(Output_dataNames(1,3))
% ylabel(Output_dataNames(1,22))
% legend(Output_dataNames(1,22))

%% dependent on shape of det cell size so it shows nothign
% figure("Name","Minimum Channel Width @ P0=5atm, T0=313.15K")
% plt=loglog(Output(:,3),Output(:,23),LineStyle="--",Marker="o");
% grid on
% 
% xlabel(Output_dataNames(1,3))
% ylabel(Output_dataNames(1,23))
% legend(Output_dataNames(1,23))

%% dependent on shape of det cell size so it shows nothign
% figure("Name","OD: P0=5atm, T0=313.15K")
% plt=loglog(Output(:,3),Output(:,24),LineStyle="--",Marker="o");
% grid on
% 
% xlabel(Output_dataNames(1,3))
% ylabel(Output_dataNames(1,24))
% legend(Output_dataNames(1,24))






















