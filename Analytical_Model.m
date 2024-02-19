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

P1 = 10000; % [Pa]
T1 = 293;% [K]
eq=1.0;
mech = 'Burke2012.yaml';
% mech = 'Hong2011.yaml';
% mech = 'sandiego20161214_H2only.yaml'; % this is fucked

% New mechanism files on SDToolbox site, Burke2012 most recent one created.
% Need to review all three papers, discuss what they do/whats new and conclude why to pick one vs the rest. 
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
vN_Point = vN_State(P1, T1, FAR, mech, gas1);

%% Calculating CJ State
CJ_Point = CJ_State(P1, T1, FAR, mech, gas1);

%% Calculating ZND Detonation Structure
Detonation_Structure = ZND_Structure_Shak(P1, T1, FAR, mech, gas1); %we only pull cell size from this

%% Geometry Definition
% Equations taken from: 
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "A Theoretical Review of Rotating Detonation Engines" 
% Ian J Shaw et al.

% Dimension are in millimeters, geometry calcs from hypergolic report
cell_gav=Detonation_Structure(1,22);
Minimum_Channel_OD = 40*cell_gav*1000;
Minimum_Channel_Width = 2.4*cell_gav*1000;
Minimum_Chamber_Length = 24*cell_gav*1000;
Minimum_Channel_ID = (Minimum_Channel_OD - 2*Minimum_Channel_Width);

disp([' '])
disp(['................................................................']);
disp(['Geometry Definition'])

disp([' '])
% disp(['Cell Size: ', num2str(Detonation_Structure(1,22)),' (m)']);
disp(['Minimum Channel OD (lambda): ', num2str(Minimum_Channel_OD),' (mm)']);
disp(['Minimum Channel ID: ', num2str(Minimum_Channel_ID),' (mm)']);
disp(['Minimum Channel Width (delta): ', num2str(Minimum_Channel_Width),' (mm)']);
disp(['Minimum Channel Length (L): ', num2str(Minimum_Chamber_Length),' (mm)']);
% disp(['TFinal (ZND): ', num2str(Detonation_Structure(1,14)),' (K)']);

% Using geometry calcs from big red
% Ian J Shaw et al., “A Theoretical Review of Rotating Detonation Engines”, doi: 10.5772.
big_red_minimumFillHeight=(12-5)*cell_gav*1000;
big_red_minD=28*cell_gav*1000;
big_red_min_delta=0.2*big_red_minimumFillHeight;
big_red_minLength=2*big_red_minimumFillHeight;
disp([' '])
% disp(['................................................................']);
disp([' '])
disp(['Geometry Per Big Red Rules of Thumb'])
disp(['Minimum Fill height: ', num2str(big_red_minimumFillHeight),' (mm)'])
disp(['Minimum Diameter: ', num2str(big_red_minD),' (mm)'])
disp(['Minimum Delta: ', num2str(big_red_min_delta),' (mm)'])
disp(['Minimum Length: ', num2str(big_red_minLength),' (mm)'])


%% M_dot calculation
% Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "Rotating Detonation Wave Stability" [Piotr Wolanski]
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]

Fill_Height = (12-5)*cell_gav; %This is the max critical fill height case

% LOGAN HARD CODING GEOMETRY SELECTED TO GET MDOT
Minimum_Channel_OD=55.11;
Minimum_Channel_Width=4.70;
Minimum_Channel_ID=Minimum_Channel_OD-2*Minimum_Channel_Width;

% Fill_Height = (12)*cell_gav;

Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2)-((Minimum_Channel_ID/1000).^2))*Fill_Height;

m_dot_tot_exit = Fill_Height*(Minimum_Channel_Width/1000)*CJ_Point(1,4)*CJ_Point(1,1); %density of combustion products, and cj speed
% [J. Shepherd, J. Kasahara]...can't find the equation anymore to verify the density and speed conditions

C=pi()*(Minimum_Channel_OD/1000);
m_dot_intoengine=(Fill_Volume*R1*CJ_Point(1,1)/C); 

m_dot_set = 0.32;

% m_dot_tot_exit = (Minimum_Channel_Width/1000)*C*CJ_Point(1,4)*CJ_Point(1,1);
% typical rho*v*A equation at the exit surface, though CJ speed is not the axial speed of expanded gas

% m_dot_H2 = m_dot_tot*H2_Percent;
% m_dot_O2 = m_dot_tot*O2_Percent;

disp([' '])
disp(['................................................................']);
disp(['Mass Flow'])

disp([' '])
disp(['Fill Height: ', num2str(Fill_Height),' (m)']);
disp(['Fill Vol: ', num2str(Fill_Volume),' (m^3)']);
disp(['m_dot Total (Shepherd, Kasahara): ', num2str(m_dot_tot_exit),' (kg/s)']);
disp(['m_dot Total (Us): ', num2str(m_dot_intoengine),' (kg/s)']);

%% Wave Number Calculation

Critical_Fill_Height = 7*Detonation_Structure(1,22);

Minimum_Channel_ID = Minimum_Channel_OD - Minimum_Channel_Width;
Mean_Channel_Diam = ((Minimum_Channel_OD - Minimum_Channel_ID)/2)+(Minimum_Channel_ID);

Critical_Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2) - ((Minimum_Channel_ID/1000).^2))*Critical_Fill_Height;

Fill_Time = pi*(Mean_Channel_Diam/1000)/CJ_Point(1,1);

Volumetric_Mixture_Supply = m_dot_tot_exit*CJ_Point(1,5);
t_mf = Critical_Fill_Volume/Volumetric_Mixture_Supply;

Wave_Number = Fill_Time/t_mf;



disp([' '])
disp(['................................................................']);
disp(['Wave Number'])

disp([' '])
disp(['Critical Fill Height: ', num2str(Critical_Fill_Height),' (m)']);
disp(['Critical Fill Volume: ', num2str(Critical_Fill_Volume),' (m^3)']);
disp(['Fill Time: ', num2str(Fill_Time),' (s)']);
disp(['Volumetric Mixture Supply: ', num2str(Volumetric_Mixture_Supply),' (m^3/s)']);
disp(['Wave Number: ', num2str(Wave_Number),' ']);

%% Thrust Calculation (Check)
% Equations taken from:
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]
% - "Techniques for the Estimation of Heats of Explosion (HEX) Using Thermochemical Codes" [Robert A. Fifer, Jeffery B. Morris]
% - "Numerical Solution Methods for Shock and Detonation Jump Conditions" [S. Browne, J. Ziegler, and J. E. Shepherd]

V_w = CJ_Point(1,1)/(((Minimum_Channel_OD/1000)+(Minimum_Channel_ID/1000))/2); %cj speed
f = V_w/(2*pi);
T = 1/f;
P_a = 101325;

% q_hc = -241800*m_dot_H2*Detonation_Structure(1,23)/0.00202; %heat released based on H2 heat of combustion and the amount of H2 being used
% q_HEX = ((241800-67.63*0.5*2)/2.01568)*1000; % [J/kg] final units, numerator units [J/mol] (heat of combustion and other junk), denom units [g/mol] (of compound)
% q_HEX = ((0.5*285830)/2.01568)*1000; % [J/kg]
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
T_e = m_dot_tot_exit*(sqrt(2*cp1*T1))*(sqrt(1 + Term_1b - Term_2b * Term_3b * Term_4b));

Tsp_e = T_e/m_dot_tot_exit;

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
%Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]

% Isp_ue = T_ue/(m_dot_tot*9.81);

Isp_e = T_e/(m_dot_intoengine*9.81);

disp([' '])
disp(['................................................................']);
disp(['Specific Impulse'])

disp([' '])
% disp(['Isp (Under Expanded): ', num2str(Isp_ue),' (s)']);
disp(['Isp (Expanded): ', num2str(Isp_e),' (s)']);
disp(['Isp (test with mdot IN): ', num2str(T_e/(m_dot_intoengine*9.81)),' (s)']);
 %% Fill Parameters
% % Equations taken from: 
% % - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% % Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% % - "Rotaing Detonation Wave Stability" [Piotr Wolanski]
% 
% Fill_Parameters = [];
% fillParameters=[];
% sz=0;
% 
% cellsize_multiple=linspace(7,17,11); % increase the third number to get finer data resolution.
% for i = cellsize_multiple
%     %%%%%%% @shak i changed the way you were adding data to the
%     %%%%%%% Fill_Parameters variable so they were divided by column; I know
%     %%%%%%% what you were going for but you accidentally had them being
%     %%%%%%% added to funky places
%     sz=size(fillParameters);
%     % Fill  Height (in m)
%     Critical_Fill_Height = i*Detonation_Structure(1,22);
% %     Fill_Parameters(end+1) = Critical_Fill_Height;
%     
%     % Fill Volume (in  m^3)
%     Minimum_Channel_ID = Minimum_Channel_OD - Minimum_Channel_Width;
%     Mean_Channel_Diam = ((Minimum_Channel_OD - Minimum_Channel_ID)/2)+(Minimum_Channel_ID);
%     
%     Critical_Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2) - ((Minimum_Channel_ID/1000).^2))*Critical_Fill_Height;
% %     Fill_Parameters(end+1) = Critical_Fill_Volume;
%     
%     % Fill Time (in s)
%     Fill_Time = pi*(Mean_Channel_Diam/1000)/CJ_Point(1,1);
% %     Fill_Parameters(end+1) = Fill_Time;
%     
%     % Detonation Wave Number
%     Volumetric_Mixture_Supply = m_dot_eq*CJ_Point(1,5);
%     t_mf = Critical_Fill_Volume/Volumetric_Mixture_Supply;
%     
%     Wave_Number = Fill_Time/t_mf;
% %     Fill_Parameters(end+1) = Wave_Number;
%     
%     fillParameters(sz(1,1)+1,:)=[Critical_Fill_Height,Critical_Fill_Volume,Fill_Time,Wave_Number]; %add data to new fill parameters variable once per iteration before the variables get cleared next iteration.
% 
% end
% 
% % Fill_Parameters = reshape(Fill_Parameters,[11,4]);
% % [row,col] = find(Fill_Parameters==max(Fill_Parameters));
% 
% disp([' '])
% disp(['................................................................']);
% disp("Fill parameters")
% % disp(Fill_Parameters)
% disp(['Based on Cell Size ', num2str(Detonation_Structure(1,22)),' (mm)']);
% % disp(fillParameters)
% 
% figure("Name","Critical Params Based on Cell Size Multiple")
% scatter(cellsize_multiple,fillParameters(:,1),DisplayName="Critical Fill Height [m]")
% hold on
% scatter(cellsize_multiple,fillParameters(:,2),DisplayName="Critical Fill Volume [m^3??]")
% scatter(cellsize_multiple,fillParameters(:,3),DisplayName="Fill Time [s]")
% xlabel("Cell size multiple (i think)")
% ylabel("Various")
% 
% yyaxis right
% scatter(cellsize_multiple,fillParameters(:,4),DisplayName="Wave Number")
% ylabel("Wave Number")
% legend
% 
% 
