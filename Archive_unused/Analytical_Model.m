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
close all force
clear
clc
disp('Analytical_Model')

P1 = 130e+3; % [Pa]
T1 = 300;% [K]
eq = 1.0;
mech = 'Burke2012.yaml';
% mech = 'h2o2.yaml';
% mech = 'Hong2011.yaml';
% mech = 'sandiego20161214_H2only.yaml'; % this is fucked

gas1 = Solution(mech);
eq=InitialState(T1,P1,eq,gas1); %%this calculates the mol ratio of hydrogen to oxygen
FAR = sprintf('H2:%d O2:%d',eq(1,1),eq(2,1));
H2_Percent = eq(1,1);
O2_Percent = eq(2,1);

% ALTERNATE -- AIR
% mech = 'airNASA9ions.yaml'; %%yaml files come from here: C:\Program Files\Cantera\data
% FAR = 'O2:0.2095 N2:0.7808 CO2:0.0004 Ar:0.0093';

%% Calculate & Print Initial State
R1 = density(gas1); %[kg/m^3]
c1_fr = soundspeed_fr(gas1); %[m/s]
cp1 = cp_mass(gas1); % specific heat capacity of 1kg of the mixture [J/kg* K]
w1 = meanMolecularWeight(gas1); %mean molecular weight [kg/kmol]
gamma1_fr =  c1_fr*c1_fr*R1/P1; %(v^2)*rho/P %ratio of heat capacities (cp/cv)
H1 = enthalpy_mass(gas1);
alt_gamma_1 = cp_mass(gas1)/cv_mass(gas1); % gamma is the ratio of heat capacities; we have no clue 
    %how the original gamma1_fr calculates the right answer; units dont make sense but alt_gamma_1 validates the solution

disp([' ']);
disp(['................................................................']);
disp(['Initial state']);

disp([' ']);
disp(['   Pressure: ',num2str(P1),' (Pa)']);
disp(['   Temperature: ',num2str(T1),' (K)']);
disp(['   Density: ',num2str(R1),' (kg/m3)']);
disp(['   Enthalpy: ',num2str(H1),' (kJ/kg)']);
disp(['   c1 (frozen): ',num2str(c1_fr),' (m/s) -- speed of sound in propellant']); %
disp(['   gamma1 (frozen): ',num2str(gamma1_fr),'(unitless)']); %
disp(['   Cp1 (frozen): ',num2str(cp1),' (J/kg*K) ']); %
disp(['   Mean Molecular Weight: ',num2str(w1),'(kg/kmol) ']); 

%% Calculating von Neumann Point
vN_Point = vN_State(P1, T1, FAR, mech, gas1,1);

%% Calculating CJ State
CJ_Point = CJ_State(P1, T1, FAR, mech, gas1,1);

%% Calculating ZND Detonation Structure
Detonation_Structure = ZND_Structure_Shak(P1, T1, FAR, mech, gas1); %we only pull cell size from this

cell_sean = ((1.6e-3*101325)/P1);
disp(['Cell Size (Sean): ', num2str(cell_sean),' (m)']);
%% Geometry Definition
% Equations taken from: 
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "A Theoretical Review of Rotating Detonation Engines" 
% Ian J Shaw et al.

% Dimension are in millimeters, geometry calcs from hypergolic report
cell_west = 29*Detonation_Structure(1,18);
Minimum_Channel_OD = 40*cell_sean;
Minimum_Channel_Width = 2.4*cell_sean;
Minimum_Chamber_Length = 24*cell_sean;
Minimum_Channel_ID = (Minimum_Channel_OD - (2*Minimum_Channel_Width));

disp([' '])
disp(['................................................................']);
disp(['Geometry Definition'])

disp([' '])
disp(['Minimum Channel OD (lambda): ', num2str(Minimum_Channel_OD),' (mm)']);
disp(['Minimum Channel ID: ', num2str(Minimum_Channel_ID),' (mm)']);
disp(['Minimum Channel Width (delta): ', num2str(Minimum_Channel_Width),' (mm)']);
disp(['Minimum Channel Length (L): ', num2str(Minimum_Chamber_Length),' (mm)']);

% % Using geometry calcs from big red
% % Ian J Shaw et al., “A Theoretical Review of Rotating Detonation Engines”, doi: 10.5772.
% big_red_minimumFillHeight=(12-5)*cell_west*1000;
% big_red_minD=28*cell_west*1000;
% big_red_min_delta=0.2*big_red_minimumFillHeight;
% big_red_minLength=2*big_red_minimumFillHeight;
% disp([' '])
% % disp(['................................................................']);
% disp([' '])
% disp(['Geometry Per Big Red Rules of Thumb'])
% disp(['Minimum Fill height: ', num2str(big_red_minimumFillHeight),' (mm)'])
% disp(['Minimum Diameter: ', num2str(big_red_minD),' (mm)'])
% disp(['Minimum Delta: ', num2str(big_red_min_delta),' (mm)'])
% disp(['Minimum Length: ', num2str(big_red_minLength),' (mm)'])


%% M_dot calculation
% Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "Rotating Detonation Wave Stability" [Piotr Wolanski]
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]

Cl = 12;
Fill_Height = (Cl)*cell_sean; %This is the max critical fill height case

Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2)-((Minimum_Channel_ID/1000).^2))*Fill_Height;

m_dot_P_history = Fill_Height*(Minimum_Channel_Width)*density(gas1)*CJ_Point(1,1); %density of combustion products, and cj speed
% [J. Shepherd, J. Kasahara]

C = pi()*(Minimum_Channel_OD/1000);
m_dot_intoengine=(Fill_Volume*CJ_Point(1,4)*CJ_Point(1,1)/C); 

% m_dot_H2 = m_dot_tot*H2_Percent;
% m_dot_O2 = m_dot_tot*O2_Percent;

disp([' '])
disp(['................................................................']);
disp(['Mass Flow'])

disp([' '])
disp(['Fill Height: ', num2str(Fill_Height),' (m)']);
disp(['Fill Vol: ', num2str(Fill_Volume),' (m^3)']);
disp(['m_dot Total (Shepherd, Kasahara): ', num2str(m_dot_P_history),' (kg/s)']);
disp(['m_dot Total (Us): ', num2str(m_dot_intoengine),' (kg/s)']);

%% Wave Number Calculation
% Equations taken from:
% - "Rotaing Detonation Wave Stability" [Piotr Wolanski]
% - "Small-size rotating detonation engine: scaling and minimum mass flow
% rate" [Sean Connolly-Boutin et al]

R_sp = 8.314462618;

Critical_Fill_Height = Fill_Height;

Minimum_Channel_ID = Minimum_Channel_OD - Minimum_Channel_Width;
Mean_Channel_Diam = ((Minimum_Channel_OD - Minimum_Channel_ID)/2)+(Minimum_Channel_ID);

Critical_Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2) - ((Minimum_Channel_ID/1000).^2))*Critical_Fill_Height;

Fill_Time = pi*(Mean_Channel_Diam/1000)/CJ_Point(1,1);

Volumetric_Mixture_Supply = m_dot_P_history*CJ_Point(1,5);
t_mf = Critical_Fill_Volume/Volumetric_Mixture_Supply;

Wave_Number = Fill_Time/t_mf;

Wave_Number_Sean = (m_dot_P_history*R_sp*T1)/(Cl*0.0016*101325*CJ_Point(1,1)*Critical_Fill_Height);

Fill_Time_Sean = (pi*(Mean_Channel_Diam))/(CJ_Point(1,1)*Wave_Number_Sean);

disp([' '])
disp(['................................................................']);
disp(['Wave Number'])

disp([' '])
disp(['Critical Fill Height: ', num2str(Critical_Fill_Height),' (m)']);
disp(['Critical Fill Volume: ', num2str(Critical_Fill_Volume),' (m^3)']);
disp(['Fill Time: ', num2str(Fill_Time),' (s)']);
disp(['Volumetric Mixture Supply: ', num2str(Volumetric_Mixture_Supply),' (m^3/s)']);
disp(['Wave Number (Wolanski): ', num2str(Wave_Number),' ']);
disp([' '])
disp(['Wave Number (Sean): ', num2str(Wave_Number_Sean),' ']);
disp(['Fill Time (Sean): ', num2str(Fill_Time_Sean),' (s)']);


%% Thrust Calculation (Check)
% Equations taken from:
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]
% - "Numerical Solution Methods for Shock and Detonation Jump Conditions" [S. Browne, J. Ziegler, and J. E. Shepherd]

V_w = CJ_Point(1,1)/(((Minimum_Channel_OD/1000)+(Minimum_Channel_ID/1000))/2); %cj speed
f = V_w/(2*pi);
T = 1/f;
P_a = 101325;

q_h_0 = (241820/18.01528)*1000; % [J/kg] same method as used in hugoniot spread sheet, *1000 is for g->kg conversion

q_h_1 = (c1_fr^2)*((((1+gamma1_fr*(((CJ_Point(1,1)/c1_fr)^2)))^2)/(2*((CJ_Point(1,14)^2)-1)))*((CJ_Point(1,14)/gamma1_fr)^2)*(1/((CJ_Point(1,1)/c1_fr)^2))-(1/(gamma1_fr-1))-(((CJ_Point(1,1)/c1_fr)^2)/2));
    % Where: vn_Point(1,13) == M1
           % CJ_Point(1,14) == gamma2
           % CJ_Point(1,1) == CJ Speed
    % From SDToolbox manual
           
q_h = ((c1_fr^2)/(2*((CJ_Point(1,14)^2)-1)))*(((CJ_Point(1,1)/c1_fr)-(1/(CJ_Point(1,1)/c1_fr)))^2);
    % From SDToolbox manual

% All three methods uncommented above of Q now currently align (enough)

T_t = 1350; % [N]

Term_1b = (q_h)/(cp1*T1);
Term_2b = (P_a/P1).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
Term_3b = (P1/CJ_Point(1,2)).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
Term_4b = (CJ_Point(1,3)/T1);
T_e = m_dot_P_history*(sqrt(2*cp1*T1))*(sqrt(1 + Term_1b - Term_2b * Term_3b * Term_4b));

Tsp_e = T_e/m_dot_P_history;

disp([' '])
disp(['................................................................']);
disp(['Thrust'])

disp(['Angular Velocity: ', num2str(V_w),' (rad/s)']);
disp(['Wave Frequency: ', num2str(f),' (Hz)']);
disp([' '])
disp(['Target Thrust: ', num2str(T_t),' (N)']);
% disp(['Actual Thrust (Under Expanded): ', num2str(T_ue), ' (N)']);
disp(['Thrust (Expanded): ', num2str(T_e), ' (N)']);
% disp(['Specific Thrust (Under Expanded): ', num2str(Tsp_ue), ' (N/kg/s)']);
disp(['Specific Thrust (Expanded): ', num2str(Tsp_e), ' (N/kg/s)']);

%% Isp Calculation
% Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]

Isp_e = T_e/(m_dot_P_history*9.81);

disp([' '])
disp(['................................................................']);
disp(['Specific Impulse'])

disp([' '])
disp(['Isp (Expanded): ', num2str(Isp_e),' (s)']);
disp(['Isp: ', num2str(T_e/(m_dot_P_history*9.81)),' (s)']);
