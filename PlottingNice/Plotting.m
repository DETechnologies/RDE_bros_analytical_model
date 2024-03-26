close all force 
clear 
clc

load("AnalyticalModel_calculatorOutput_vPressure_March7_fixed.mat")

%thrust plot ---------------------------------------------------------------
figure("Name","pressure vs thrust")
plot(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Thrust O/P (N)"});
hold on
plot([24e+3;880e+3],[1350;1350],Color="black"); % plot horizontal line at 1350N
xlim([24e+3 880e+3])
xlabel("Initial Pressure P_0 [Pa]")
ylabel("Thrust [N]")

plot([101.325e+3;101.325e+3],[0;2500],Color="black",LineStyle="--"); %this plots a vertical line showing where atmostpheric pressure is.
yyaxis right
ylabel("Detonation Cell Size \lambda [mm]")

loglog(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Cell Size value (m)"},Color="black",LineStyle=":",LineWidth=2);
legend("Thrust","Thrust Objective (1350N)","Atmospheric Pressure (101.325kPa)","SeanCB Cell Size")

%mdot plot ---------------------------------------------------------------
figure("Name","pressure vs mdot")
plot(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"mDot (kg/s)"});
xlabel("Initial Pressure P_0 [Pa]")
ylabel("mass flow rate")

%UCJ and cell size plot
figure("Name","pressure vs mdot")
plot(Outputs{:,"CJ Speed (m/s)"},Outputs{:,"Cell Size value (m)"});
ylabel("Cell Size")
xlabel("CJ Speed")

%plot all cell sizes
detonationDatabase=readtable("DetonationDatabase_data.xlsx");
detonationDatabase.Pressure_kpa_ = detonationDatabase.Pressure_kpa_*1000;
detonationDatabase.CellSize_mm_ = detonationDatabase.CellSize_mm_/1000;
figure("Name","pressure vs mdot")
loglog(Outputs{:,"I/P Pressure (Pa)"},[Outputs{:,"SeanCB cell Size (m)"},Outputs{:,"Westbrook Cell Size (m)"},Outputs{:,"Gav Cell Size (m)"},Outputs{:,"NG Cell Size (m)"}]);
hold on
scatter(table2array(detonationDatabase(:,1)),table2array(detonationDatabase(:,2)),'x','black');
% legend("seancb","west","gav","ng")
grid on
set(gca,'xscale','log','yscale','log')
ylabel("Cell Size")
xlabel("IP Pressure")

legend("CalTech Detonation Database","Westbrook","Garikov","Ng et al.","Sean CB")

