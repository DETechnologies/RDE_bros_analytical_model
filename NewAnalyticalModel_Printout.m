% The new analytical model
% DETechnologies
% Logan and Shak - 2023/2024

close all force 
clear 
clc

%% init params
P1 = 130e+3; % [Pa]
T1 = 300;% [K]
eq = 1.0;
mechFiles={'Burke2012.yaml';'h2o2.yaml';'Hong2011.yaml'};
mech = mechFiles{2}; %1 means it uses the first one (burke) (first row (row matrix))

%% Call function that does the actual maths
[gas1,vN_Point,CJ,ZND,CellSizePredictions,Misc,GeometryPredictor] = NewAnalyticalModel(P1,T1,eq,mech,0,4,1,1);
% % % % % % % % -5 < Bykovskii_adder < 5
% % % % % % % % CellCorr2Use one of these NUMBERS {'Gavrikov','Westbrook','Ng','SeanCB'} = 1,2,3,4
% % % % % % % % Geometry rule; 0 = nair, 1 = bykovskii.
% % % % % % % % Print things 1/0

%% Printout Section - for shak

%printing initial state things here.
fprintf('\n\nInitial State')
fprintf('\nPressure: %d [Pa]',P1); fprintf('\nTemperature: %d [K]',T1); fprintf('\nEquivalence Ratio: %d',eq)
fprintf('\nDensity (initial): %d [kg/m3]',density(gas1)); fprintf('\nEnthalpy: %d [kJ/kg]',density(gas1));
fprintf('\nSpeed of Sound (premixed propellant): %d',soundspeed_fr(gas1));

%printing out geometry table
fprintf('\n\nGeometry')
GeometryPredictor%dont hide output (no semicolon!!)

%Printing thrust
fprintf('\n\nPerformance Indicators')
fprintf('\nThrust: %d [N]',Misc{1,"Thrust"})
fprintf('\nSpecific Impulse: %d [s^-1]',Misc{1,"ISP"})
fprintf('\nFill Time: %d [s]',Misc{1,"Fill_Time_Sean"})
fprintf('\nWave Number: %d',Misc{1,"Wave_Number_Sean"})
fprintf('\nThrust Goal: %d [N]', Misc{1,"Thrust_Goal"})

%Printing mass flow rate
fprintf('\n\nMass Flow')
fprintf('\nThrust Based Mass Flow: %d [kg/s]', Misc{1,"m_dot_T"})
fprintf('\nChamber Volume Mass Flow: %d [kg/s]', Misc{1,"m_dot_V"})
fprintf('\nPropellant Fill Height to hit 1350N: %d [m]',Misc{1,"m_dot_T"}/(1*5e-3*density(gas1)*CJ(1)))
