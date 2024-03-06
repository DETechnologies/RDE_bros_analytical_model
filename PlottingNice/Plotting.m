close all force 
clear 
clc

load("AnalyticalModel_calculatorOutput_vPressure_March5.mat")

%thrust plot ---------------------------------------------------------------
figure("Name","pressure vs thrust")
plot(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Thrust O/P"});
hold on
plot([24e+3;880e+3],[1350;1350],Color="black"); % plot horizontal line at 1350N
xlim([24e+3 880e+3])
xlabel("Initial Pressure P_0 [Pa]")
ylabel("Thrust [N]")

plot([101.325e+3;101.325e+3],[1000;2500],Color="black",LineStyle="--"); %this plots a vertical line showing where atmostpheric pressure is.
yyaxis right
ylabel("Detonation Cell Size \lambda [mm]")

loglog(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Cell Size value"},Color="black",LineStyle=":",LineWidth=2);
legend("Thrust","Thrust Objective (1350N)","Atmospheric Pressure (101.325kPa)","SeanCB Cell Size")



%mdot plot ---------------------------------------------------------------
figure("Name","pressure vs mdot")
plot(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"mDot"});
xlabel("Initial Pressure P_0 [Pa]")
ylabel("mass flow rate")



%UCJ and cell size plot
figure("Name","pressure vs mdot")
plot(Outputs{:,"CJ Speed (m/s)"},Outputs{:,"Cell Size value"});
ylabel("Cell Size")
xlabel("CJ Speed")

