% ################################################################################
% Theory, numerical methods and applications are described in the following report:
% 
%     Numerical Solution Methods for Shock and Detonation Jump Conditions, S.
%     Browne, J. Ziegler, and J. E. Shepherd, GALCIT Report FM2006.006 - R3,
%     California Institute of Technology Revised September, 2018.
% 
% Please cite this report and the website if you use these routines. 
% 
% Please refer to LICENCE.txt or the above report for copyright and disclaimers.
% 
% http://shepherd.caltech.edu/EDL/PublicResources/sdt/
% 
% ################################################################################
%% Set Input Parameters
clear;
clc;
close all;
disp('Analytical_Model')

P1 = 2473.+3; % [Pa]
T1 = 323.15; % [K]
eq=1.15;
mech = 'h2o2.yaml'; %%yaml files come from here: C:\Program Files\Cantera\data
gas1 = Solution(mech);
q=InitialState(T1,P1,eq,gas1); %%this calculates the mol ratio of hydrogen to oxygen
% set(gas1,'Temperature',T1,'Pressure',P1,'MoleFractions',q);

%% Calculate & Print Initial State
R1 = density(gas1); %[kg/m^3]
c1_fr = soundspeed_fr(gas1); %[m/s]
cp1 = cp_mass(gas1); %heat capacity of 1kg of the mixture [J]
w1 = meanMolecularWeight(gas1); %mean molecular weight [kg/kmol]
gamma1_fr =  c1_fr*c1_fr*R1/P1; %(v^2)*rho/P %ratio of heat capacities (cp/cv)
alt_gamma_1 = cp_mass(gas1)/cv_mass(gas1); % gamma is the ratio of heat capacities; we have no clue 
    %how the original gamma1_fr calculates the right answer; units dont make sense but alt_gamma_1 validates the solution

disp([' ']);
disp(['................................................................']);
disp(['Initial state']);
disp([' ']);
disp(['   Pressure: ',num2str(P1),' (Pa)']);
disp(['   Temperature: ',num2str(T1),' (K)']);
disp(['   Density: ',num2str(R1),' (kg/m3)']);
disp(['   c1 (frozen): ',num2str(c1_fr),' (m/s) -- speed of sound in propellant']); %
disp(['   gamma1 (frozen): ',num2str(gamma1_fr),'(unitless)']); %
disp(['   Cp1 (frozen): ',num2str(cp1),' (J/K) ']); %
disp(['   Mean Molecular Weight:  ',num2str(w1),'(kg/kmol) ']); 

%% Calculating von Neumann Point
vN_Point = vN_State_Shak(P1, T1, q, mech, gas1);

%% Calculating CJ State
CJ_Point = CJ_State_Shak(P1, T1, q, mech, gas1);

%% Calculating ZND Detonation Structure
Detonation_Structure = ZND_Structure_Shak(P1, T1, q, mech, gas1);

%% Calculation of Rocket Impulse
Impulse = Rocket_Impulse_Shak(T1, P1, q, mech);

%% M_dot calculation
% Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]

Thrust = 1000;
g = 9.81;
Isp_fr = Impulse(1,9);
Isp_eq = Impulse(1,8);
m_dot_fr = (Thrust)/(Isp_fr*g);
m_dot_eq = (Thrust)/(Isp_eq*g);

disp([' '])
disp(['................................................................']);
disp(['Mass Flow'])

disp([' '])
disp(['Set Thrust = ', num2str(Thrust),' (N)']);
disp(['M_dot (frozen)', num2str(m_dot_fr),' (kg/s)']);
disp(['M_dot (equilibrium)',  num2str(m_dot_eq),' (kg/s)']);

%% Geometry Definition
% Equations taken from: 
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]

% Dimension are in millimeters
Minimum_Channel_OD = 40*Detonation_Structure(1,22)*1000;
Minimum_Channel_Width = 2.4*Detonation_Structure(1,22)*1000;
Minimum_Chamber_Length = 24*Detonation_Structure(1,22)*1000;

disp([' '])
disp(['................................................................']);
disp(['Minimum Channel OD (lambda)', num2str(Minimum_Channel_OD),' (mm)']);
disp(['Minimum Channel Width (lambda)', num2str(Minimum_Channel_Width),' (mm)']);
disp(['Minimum Channel Length (lambda)', num2str(Minimum_Chamber_Length),' (mm)']);

%% Fill Parameters
% Equations taken from: 
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "Rotaing Detonation Wave Stability" [Piotr Wolanski]

Fill_Parameters = [];

for i = 7:1:17
% Fill  Height (in m)
Critical_Fill_Height = i*Detonation_Structure(1,22);
Fill_Parameters(end+1) = Critical_Fill_Height;

% Fill Volume (in  m^3)
Minimum_Channel_ID = Minimum_Channel_OD - Minimum_Channel_Width;
Mean_Channel_Diam = ((Minimum_Channel_OD - Minimum_Channel_ID)/2)+(Minimum_Channel_ID);

Critical_Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2) - ((Minimum_Channel_ID/1000).^2))*Critical_Fill_Height;
Fill_Parameters(end+1) = Critical_Fill_Volume;

% Fill Time (in s)
Fill_Time = pi*(Mean_Channel_Diam/1000)/CJ_Point(1,1);
Fill_Parameters(end+1) = Fill_Time;

% Detonation Wave Number
Volumetric_Mixture_Supply = m_dot_eq*CJ_Point(1,5);
t_mf = Critical_Fill_Volume/Volumetric_Mixture_Supply;

Wave_Number = Fill_Time/t_mf;
Fill_Parameters(end+1) = Wave_Number;

end

Fill_Parameters = reshape(Fill_Parameters,[11,4]);
[row,col] = find(Fill_Parameters==max(Fill_Parameters));

disp([' '])
disp(['................................................................']);

