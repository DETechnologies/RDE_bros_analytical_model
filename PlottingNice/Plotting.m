close all force 
clear 
clc

% load("AnalyticalModel_calculatorOutput_vPressure_March7_fixed.mat")
load("combined.mat")

%thrust plot ---------------------------------------------------------------
figure("Name","pressure vs thrust")
plot(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Thrust O/P (N)"});
hold on
plot([0e+3;880e+3],[1350;1350],Color="black"); % plot horizontal line at 1350N
xlim([10e+3 700e+3])
ylim([0,2500])
xlabel("Initial Pressure P_0 [Pa]")
ylabel("Thrust [N]")

plot([101.325e+3;101.325e+3],[0;3500],Color="black",LineStyle="--"); %this plots a vertical line showing where atmostpheric pressure is.
yyaxis right
ylabel("Detonation Cell Size \lambda [mm]")

loglog(Outputs{:,"I/P Pressure (Pa)"},Outputs{:,"Cell Size value (m)"},Color="black",LineStyle=":",LineWidth=2);
legend("Thrust","Thrust Objective (1350N)","Atmospheric Pressure (101.325kPa)","Connolly-Boutin Cell Size")

%plot all cell sizes -------------------------------------------------------
detonationDatabase=readtable("DetonationDatabase_data.xlsx");
detonationDatabase.Pressure_kpa_ = detonationDatabase.Pressure_kpa_*1000;
detonationDatabase.CellSize_mm_ = detonationDatabase.CellSize_mm_/1000;
figure("Name","pressure vs mdot")
loglog(Outputs{:,"I/P Pressure (Pa)"},[Outputs{:,"SeanCB cell Size (m)"},Outputs{:,"Westbrook Cell Size (m)"},Outputs{:,"Gav Cell Size (m)"},Outputs{:,"NG Cell Size (m)"}]);
hold on
scatter(table2array(detonationDatabase(:,1)),table2array(detonationDatabase(:,2)),'x','black');

%plot line of best fit
% [exp_tr,gof_tr] = fit(table2array(detonationDatabase(:,1)),table2array(detonationDatabase(:,2)),"exp2",Algorithm="Levenberg-Marquardt");
% plot(exp_tr)

% grid on
set(gca,'xscale','log','yscale','log')
ylabel("Cell Size [m]")
xlabel("IP Pressure [Pa]")
xlim([4e+3,1.2e+6])

legend("Connolly-Boutin","Westbrook","Gavrikov","Ng et al.","CalTech Detonation Database")

% plot geometry functions
Bykov_hstar_lo=12-5;
Bykov_hstar_hi=12+5;
Bykov_D=28;
BykovDelta=1/5;
Bykov_L_lo=2*(12-5);
Bykov_L_hi=2*(12+5);
Bykov=table(Bykov_hstar_lo,Bykov_hstar_hi,Bykov_D,BykovDelta,Bykov_L_lo,Bykov_L_hi);

nair_hstar_lo=12-5;
nair_hstar_hi=12+5;
nair_D=40;
nairDelta=2.4;
nair_L=24;
Nair=table(nair_hstar_lo,nair_hstar_hi,nair_D,nairDelta,nair_L);


% now multiply cell size correlations by each.
GeometryData=table(Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_hstar_lo, Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_hstar_hi, Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_D, Outputs.("SeanCB cell Size (m)")*Bykov.BykovDelta,Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_L_lo,Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_L_hi, ...
                   Outputs.("SeanCB cell Size (m)")*Nair.nair_hstar_lo, Outputs.("SeanCB cell Size (m)")*Nair.nair_hstar_hi, Outputs.("SeanCB cell Size (m)")*Nair.nair_D, Outputs.("SeanCB cell Size (m)")*Nair.nairDelta,Outputs.("SeanCB cell Size (m)")*Nair.nair_L, ...
                   Outputs.("Westbrook Cell Size (m)")*Bykov.Bykov_hstar_lo, Outputs.("Westbrook Cell Size (m)")*Bykov.Bykov_hstar_hi, Outputs.("Westbrook Cell Size (m)")*Bykov.Bykov_D, Outputs.("Westbrook Cell Size (m)")*Bykov.BykovDelta,Outputs.("Westbrook Cell Size (m)")*Bykov.Bykov_L_lo,Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_L_hi, ...
                   Outputs.("Westbrook Cell Size (m)")*Nair.nair_hstar_lo, Outputs.("Westbrook Cell Size (m)")*Nair.nair_hstar_hi, Outputs.("Westbrook Cell Size (m)")*Nair.nair_D, Outputs.("Westbrook Cell Size (m)")*Nair.nairDelta,Outputs.("Westbrook Cell Size (m)")*Nair.nair_L, ...
                   Outputs.("NG Cell Size (m)")*Bykov.Bykov_hstar_lo, Outputs.("NG Cell Size (m)")*Bykov.Bykov_hstar_hi, Outputs.("NG Cell Size (m)")*Bykov.Bykov_D, Outputs.("NG Cell Size (m)")*Bykov.BykovDelta,Outputs.("NG Cell Size (m)")*Bykov.Bykov_L_lo,Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_L_hi, ...
                   Outputs.("NG Cell Size (m)")*Nair.nair_hstar_lo, Outputs.("NG Cell Size (m)")*Nair.nair_hstar_hi, Outputs.("NG Cell Size (m)")*Nair.nair_D, Outputs.("NG Cell Size (m)")*Nair.nairDelta,Outputs.("NG Cell Size (m)")*Nair.nair_L, ...
                   Outputs.("Gav Cell Size (m)")*Bykov.Bykov_hstar_lo, Outputs.("Gav Cell Size (m)")*Bykov.Bykov_hstar_hi, Outputs.("Gav Cell Size (m)")*Bykov.Bykov_D, Outputs.("Gav Cell Size (m)")*Bykov.BykovDelta,Outputs.("Gav Cell Size (m)")*Bykov.Bykov_L_lo,Outputs.("SeanCB cell Size (m)")*Bykov.Bykov_L_hi, ...
                   Outputs.("Gav Cell Size (m)")*Nair.nair_hstar_lo, Outputs.("Gav Cell Size (m)")*Nair.nair_hstar_hi, Outputs.("Gav Cell Size (m)")*Nair.nair_D, Outputs.("Gav Cell Size (m)")*Nair.nairDelta,Outputs.("Gav Cell Size (m)")*Nair.nair_L,...
    'VariableNames',{'SeanCB_bykov_hstar_lo','SeanCB_bykov_hstar_hi','SeanCB_bykov_D','SeanCB_bykov_delta','SeanCB_bykov_L_lo','SeanCB_bykov_L_hi',...
                    'SeanCB_nair_hstar_lo','SeanCB_nair_hstar_hi','SeanCB_nair_D','SeanCB_nair_delta','SeanCB_nair_L',...
                    'Westbrook_bykov_hstar_lo','Westbrook_bykov_hstar_hi','Westbrook_bykov_D','Westbrook_bykov_delta','Westbrook_bykov_L_lo','Westbrook_bykov_L_hi',...
                    'Westbrook_nair_hstar_lo','Westbrook_nair_hstar_hi','Westbrook_nair_D','Westbrook_nair_delta','Westbrook_nair_L',...
                    'NgetAl_bykov_hstar_lo','NgetAl_bykov_hstar_hi','NgetAl_bykov_D','NgetAl_bykov_delta','NgetAl_bykov_L_lo','NgetAl_bykov_L_hi',...
                    'NgetAl_nair_hstar_lo','NgetAl_nair_hstar_hi','NgetAl_nair_D','NgetAl_nair_delta','NgetAl_nair_L',...
                    'Gav_bykov_hstar_lo','Gav_bykov_hstar_hi','Gav_bykov_D','Gav_bykov_delta','Gav_bykov_L_lo','Gav_bykov_L_hi',...
                    'Gav_nair_hstar_lo','Gav_nair_hstar_hi','Gav_nair_D','Gav_nair_delta','Gav_nair_L'});


% take this monster and plot each correlation separately.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure("Name","Critical Fill Height")

hold on
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_hstar_lo,Color=[0, 0.4470, 0.7410])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_hstar_hi,Color=[0, 0.4470, 0.7410],LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_hstar_lo,Color=[0.6350, 0.0780, 0.1840]	)
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_hstar_hi,Color=[0.6350, 0.0780, 0.1840]	,LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_hstar_lo,Color=[0.9290, 0.6940, 0.1250])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_hstar_hi,Color=[0.9290, 0.6940, 0.1250],LineStyle="--")
 
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_hstar_lo,Color=[0.4940, 0.1840, 0.5560])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_hstar_hi,Color=[0.4940, 0.1840, 0.5560],LineStyle="--")

plot([101.325e+3 101.325e+3],[0 6],Color='black',LineStyle="--")

xlabel('Initial Pressure [Pa]')
ylabel('Critical Fill Height [m]')
title('Critical Fill Height, h*')

xlim([100e+3 200e+3])
ylim([0 0.05])

legend("Connolly-Boutin (Bykovskii) Lo","Connolly-Boutin (Bykovskii) Hi",...
    "Westbrook (Bykovskii) Lo","Westbrook (Bykovskii) Hi",...
    "Gavrikov (Bykovskii) Lo","Gavrikov (Bykovskii) Hi",...
    "Ng et al. (Bykovskii) Hi","Ng et al. (Bykovskii) Lo","P_{atm}")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(Name="Outer Diameter")

hold on

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_D,Color=[0, 0.4470, 0.7410])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_nair_D,Color=[0, 0.4470, 0.7410],LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_D,Color=[0.6350, 0.0780, 0.1840]	)
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_nair_D,Color=[0.6350, 0.0780, 0.1840]	,LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_D,Color=[0.9290, 0.6940, 0.1250])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_nair_D,Color=[0.9290, 0.6940, 0.1250],LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_D,Color=[0.4940, 0.1840, 0.5560])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_nair_D,Color=[0.4940, 0.1840, 0.5560],LineStyle="--")

plot([101.325e+3 101.325e+3],[0 6],Color='black',LineStyle="--")

xlabel('Initial Pressure [Pa]')
ylabel('Minimum Diameter, D [m]')
title('Minimum Diameter, D')

xlim([100e+3 200e+3])
ylim([0.015 0.12])

legend("Connolly-Boutin (Bykovskii)","Connolly-Boutin (Nair)","Westbrook (Bykovskii)","Westbrook (Nair)","Gavrikov (Bykovskii)","Gavrikov (Nair)","Ng et al. (Bykovskii)","Ng et al. (Nair)","P_{atm}")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(Name="Delta")
hold on

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_delta,Color=[0, 0.4470, 0.7410])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_nair_delta,Color=[0, 0.4470, 0.7410],LineStyle="--")
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_delta,Color=[0.6350, 0.0780, 0.1840]	)
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_nair_delta,Color=[0.6350, 0.0780, 0.1840]	,LineStyle="--")
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_delta,Color=[0.9290, 0.6940, 0.1250])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_nair_delta,Color=[0.9290, 0.6940, 0.1250],LineStyle="--")
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_delta,Color=[0.4940, 0.1840, 0.5560])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_nair_delta,Color=[0.4940, 0.1840, 0.5560],LineStyle="--")
plot([101.325e+3 101.325e+3],[0 6],Color='black',LineStyle="--")

xlabel('Initial Pressure [Pa]')
ylabel('Minimum Channel Width, δ [m]')
title('Minimum Channel Width, δ ')

xlim([100e+3 200e+3])
ylim([0 8e-3])

legend("Connolly-Boutin (Bykovskii)","Connolly-Boutin (Nair)","Westbrook (Bykovskii)","Westbrook (Nair)","Gavrikov (Bykovskii)","Gavrikov (Nair)","Ng et al. (Bykovskii)","Ng et al. (Nair)","P_{atm}")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(Name="Length")
hold on 

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_L_lo,Color=[0, 0.4470, 0.7410])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.SeanCB_bykov_L_hi,Color=[0, 0.4470, 0.7410],LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_L_lo,Color=[0.6350, 0.0780, 0.1840]	)
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Westbrook_bykov_L_hi,Color=[0.6350, 0.0780, 0.1840]	,LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_L_lo,Color=[0.9290, 0.6940, 0.1250])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.Gav_bykov_hstar_hi,Color=[0.9290, 0.6940, 0.1250],LineStyle="--")

plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_L_lo,Color=[0.4940, 0.1840, 0.5560])
plot(Outputs.("I/P Pressure (Pa)"),GeometryData.NgetAl_bykov_L_hi,Color=[0.4940, 0.1840, 0.5560],LineStyle="--")
plot([101.325e+3 101.325e+3],[0 6],Color='black',LineStyle="--")

plot([101.325e+3 101.325e+3],[0 6],Color='black',LineStyle="--")

xlim([100e+3 200e+3])
ylim([0 0.06])

xlabel('Initial Pressure [Pa]')
ylabel('Minimum Length [m]')
title('Minimum Length, L')

legend("Connolly-Boutin (Bykovskii) Lo","Connolly-Boutin (Bykovskii) Hi",...
    "Westbrook (Bykovskii) Lo","Westbrook (Bykovskii) Hi",...
    "Gavrikov (Bykovskii) Lo","Gavrikov (Bykovskii) Hi",...
    "Ng et al. (Bykovskii) Hi","Ng et al. (Bykovskii) Lo","P_{atm}")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(Name="PressureVCellSize")
plot(Outputs.("I/P Pressure (Pa)"),Outputs.("SeanCB cell Size (m)"),Color='black')
hold on
plot([101.325e+3 101.325e+3],[0,0.08],Color='black',LineStyle='--')
xlabel("Initial Pressure, P_o [Pa]")
ylabel("Connolly-Boutin Detonation Cell Size [m]")
title("Initial Pressure versus Detonation Cell Size")
xlim([90.000e+3,260e+3])
ylim([0.5e-3,0.002])
legend('Detonation Cell Size','Atmospheric Pressure')

vTemp=load('AnalyticalModelCalculator_April2_vTemp.mat');

figure(Name="TempVCellSize")
plot(vTemp.Outputs.("I/P Temperature (K)"),vTemp.Outputs.("SeanCB cell Size (m)"),Color='black')
xlabel("Initial Temperature, T_o [K]")
ylabel("Connolly-Boutin Detonation Cell Size [m]")
title("Initial Temperature versus Detonation Cell Size")
ylim([1.1e-3,1.4e-3])


figure(Name="TempVALLCellSizes")
plot(vTemp.Outputs.("I/P Temperature (K)"),[vTemp.Outputs{:,"SeanCB cell Size (m)"},vTemp.Outputs{:,"Westbrook Cell Size (m)"},vTemp.Outputs{:,"Gav Cell Size (m)"},vTemp.Outputs{:,"NG Cell Size (m)"}]);
xlabel("Initial Temperature, T_o [K]")
ylabel("Detonation Cell Size [m]")
title("Initial Temperature versus Detonation Cell Size")
legend("Connolly-Boutin","Westbrook","Gavrikov","Ng et al.")