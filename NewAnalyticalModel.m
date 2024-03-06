function [gas1,vN_Point,CJ,ZND,CellSizePredictions,Misc,GeometryPredictor] = NewAnalyticalModel(P1,T1,eq,mech,Bykovskii_adder,CellCorr2Use,GeometryRule,PrintThings)
% Bykovskii_adder: [-5:+5]
%       This drives the 12+5 Bykovskii relation. Used in both Bykovskii
%       geometry relations and Nair geometry relations (which are derived
%       from Bykovskii)
% CellCorr2Use one of these NUMBERS {'Gavrikov','Westbrook','Ng','SeanCB'}=1,2,3,4
%       This is only used in the thrust correlations - we need to know what
%       cell size param to grab, because its not really worth calculating
%       all of them. A potential future improvement is adding this into the
%       geometry calculating loop.
% Geometry rule; 0=nair, 1=bykovskii
%       This is only used in the thrust correlations - we need to know what
%       cell size param to grab, because its not really worth calculating
%       all of them. A potential future improvement is adding this into the
%       geometry calculating loop.
% Print things 1/0 
%       print output things to look at. Boolean.

troubleshooting=0%dont hide this line o/p --> failsafe
if troubleshooting
    %Need to make this a matlab script, rather than fcn (comment top line)
    %this inits some params to get some values
    clear
    P1 = 180e+3; % [Pa]
    T1 = 300;% [K]
    eq = 1.0;
    mech = 'Burke2012.yaml';
    CellCorr2Use=4
    GeometryRule=1
    PrintThings=1
end 

% Theory, numerical methods and applications are described in the following report:
% 
%     Numerical Solution Methods for Shock and Detonation Jump Conditions, S.
%     Browne, J. Ziegler, and J. E. Shepherd, GALCIT Report FM2006.006 - R3,
%     California Institute of Technology Revised September, 2018.
%% Initial State (pre-combustion) - State 1 (one)
CellCorr2Use=CellCorr2Use*2;
gas1 = Solution(mech);
eq=InitialState(T1,P1,eq,gas1); % mol ratio of fuel:ox
FAR = sprintf('H2:%d O2:%d',eq(1,1),eq(2,1)); % Fuel:Oxidizer mol ratio in string
P_a=101.325e+3; %P_atm in pa
%% Calculating combustion parameters to be used later
% Calculating VN Point
vN_Point = vN_State(P1, T1, FAR, mech, gas1,PrintThings);
% Calculating CJ
CJ = CJ_State(P1, T1, FAR, mech, gas1,PrintThings);
% Calculating ZND
ZND = ZND_Structure(P1, T1, FAR, mech, gas1,PrintThings); % all cell sizes in [m]
% Sean Conolly-Boutin Additions
SeanCB_CellSize = ((1.6*101.325e+3)/P1)/1000; % [m]
CellSizePredictions=table(ZND(22),ZND(18),ZND(21),SeanCB_CellSize,'VariableNames',{'Gavrikov','Westbrook','Ng','SeanCB'});
%% Geometry
% Methods:
% F. A. Bykovskii, S. A. Zhdan, and E. F. Vedernikov, "Continuous spin detonations,” Journal of propulsion and power, vol. 22, no. 6, pp. 1204–1216, 2006
% A. P. Nair, A. R. Keller, N. Q. Minesi, D. I. Pineda, and R. M. Spearrin, "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor," Proceedings of the Combustion Institute, vol. 39, no. 3, pp. 2757–2765, 2023.
if ~ exist("Bykovskii_adder",'var')
    Bykovskii_adder=0; % this is the 12+5 thing
end
for i=1:size(CellSizePredictions,2) 
    %each loop adds two rows.
    NewR=table([{sprintf('Bykovskii (12+%f)',Bykovskii_adder)};{sprintf('Nair (12+%f)',Bykovskii_adder)}],... %Bykovskii/Nair rows
        [CellSizePredictions.Properties.VariableNames(1,i);CellSizePredictions.Properties.VariableNames(1,i)],... % adds cell size predictor name (westbrook/SeanCB/Gavrikov
        [CellSizePredictions{1,i};CellSizePredictions{1,i}],... % cell size in both cols (according to loop)
        [(12+Bykovskii_adder)*CellSizePredictions{1,i};(12+Bykovskii_adder)*CellSizePredictions{1,i}],... % Min fill height
        [28*CellSizePredictions{1,i};40*CellSizePredictions{1,i}],... % min OD
        [((12+Bykovskii_adder)/5)*CellSizePredictions{1,i};2.4*CellSizePredictions{1,i}],... %min delta (channel width)
        [2*(12+Bykovskii_adder)*CellSizePredictions{1,i};24*CellSizePredictions{1,i}],... % min length
        [(12+Bykovskii_adder)*CellSizePredictions{1,i}*((12+Bykovskii_adder)/5)*CellSizePredictions{1,i}*(density(gas1))*CJ(1,1);(12+Bykovskii_adder)*CellSizePredictions{1,i}*2.4*CellSizePredictions{1,i}*(density(gas1))*CJ(1,1)],...
        'VariableNames',{'GeometryCorrelations','CellSizePredictor','CellSize','MinFillHeight','MinChannelOD','MinChannelWidth','MinChannelLength','Mass Flow Rate kg/s'});  

    if i==1
        GeometryPredictor=NewR;
    else
       GeometryPredictor=[GeometryPredictor;NewR]; 
    end
end 
%% Some Sean Things
R_sp = 8.314462618;
Wave_Number_Sean = (GeometryPredictor{CellCorr2Use-GeometryRule,'Mass Flow Rate kg/s'}*R_sp*T1)/((12+Bykovskii_adder)*0.0016*101325*CJ(1,1)*GeometryPredictor{CellCorr2Use-GeometryRule,'MinFillHeight'});
Mean_Channel_Diam=GeometryPredictor{CellCorr2Use-GeometryRule,'MinChannelOD'}-GeometryPredictor{CellCorr2Use-GeometryRule,'MinChannelWidth'};
Fill_Time_Sean = (pi*(Mean_Channel_Diam))/(CJ(1,1)*Wave_Number_Sean);

%% Thrust Calcs
q_h = ((soundspeed_fr(gas1)^2)/(2*((CJ(1,14)^2)-1)))*(((CJ(1,1)/soundspeed_fr(gas1))-(1/(CJ(1,1)/soundspeed_fr(gas1))))^2);
Term_1b = (q_h)/(cp_mass(gas1)*T1);
Term_2b = (P_a/P1).^((CJ(1,14)-1)/CJ(1,14));
Term_3b = (P1/CJ(1,2)).^((CJ(1,14)-1)/CJ(1,14));
Term_4b = (CJ(1,3)/T1);
Thrust = GeometryPredictor{CellCorr2Use-GeometryRule,'Mass Flow Rate kg/s'}*(sqrt(2*cp_mass(gas1)*T1))*(sqrt(1 + Term_1b - Term_2b * Term_3b * Term_4b));

SpecThrust = Thrust/GeometryPredictor{CellCorr2Use-GeometryRule,'Mass Flow Rate kg/s'};
ISP = Thrust/(GeometryPredictor{CellCorr2Use-GeometryRule,'Mass Flow Rate kg/s'}*9.81);

Misc=table(Wave_Number_Sean,Mean_Channel_Diam,Fill_Time_Sean,Thrust,SpecThrust,ISP);