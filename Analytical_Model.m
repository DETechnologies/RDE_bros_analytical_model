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

P1 = 850000; % [Pa]
T1 =300;% [K]
eq=1;
mech = 'h2o2.yaml'; %%yaml files come from here: C:\Program Files\Cantera\data
gas1 = Solution(mech);
eq=InitialState(T1,P1,eq,gas1); %%this calculates the mol ratio of hydrogen to oxygen
FAR = sprintf('H2:%d O2:%d',eq(1,1),eq(2,1));
H2_Percent = eq(1,1);
O2_Percent = eq(2,1);

%% Calculate & Print Initial State
R1 = density(gas1); %[kg/m^3]
c1_fr = soundspeed_fr(gas1); %[m/s]
cp1 = cp_mass(gas1); %heat capacity of 1kg of the mixture [J]
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
disp(['   Enthalpy: ',num2str(H1),' (kJ)']);
disp(['   c1 (frozen): ',num2str(c1_fr),' (m/s) -- speed of sound in propellant']); %
disp(['   gamma1 (frozen): ',num2str(gamma1_fr),'(unitless)']); %
disp(['   Cp1 (frozen): ',num2str(cp1),' (J/K) ']); %
disp(['   Mean Molecular Weight: ',num2str(w1),'(kg/kmol) ']); 

%% Calculating von Neumann Point
vN_Point = vN_State(P1, T1, FAR, mech, gas1);

%% Calculating CJ State
CJ_Point = CJ_State(P1, T1, FAR, mech, gas1);

%% Calculating ZND Detonation Structure
Detonation_Structure = ZND_Structure_Shak(P1, T1, FAR, mech, gas1);

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
disp(['Geometry Definition'])

disp([' '])
disp(['Cell Size: ', num2str(Detonation_Structure(1,22)),' (m)']);
disp(['Minimum Channel OD (lambda): ', num2str(Minimum_Channel_OD),' (mm)']);
disp(['Minimum Channel Width (lambda): ', num2str(Minimum_Channel_Width),' (mm)']);
disp(['Minimum Channel Length (lambda): ', num2str(Minimum_Chamber_Length),' (mm)']);

%% M_dot calculation
% Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]
% - "Detonation cell size of liquid hypergolic propellants: Estimation from a non-premixed combustor [Anil P. Nair,
% Alex R. Keller, Nicolas Q. Minesi, Daniel I. Pineda, R. Mitchell Spearrin]
% - "Rotaing Detonation Wave Stability" [Piotr Wolanski]
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]

Fill_Height = (12+5)*Detonation_Structure(1,22); %This is the max critical fill height case

Minimum_Channel_ID = (Minimum_Channel_OD - Minimum_Channel_Width)/1000;

Fill_Volume = 0.25*pi*(((Minimum_Channel_OD/1000).^2)-((Minimum_Channel_ID/1000).^2))*Fill_Height;

m_dot_tot = Fill_Height*(Minimum_Channel_Width/1000)*CJ_Point(1,4)*CJ_Point(1,1);
% [J. Shepherd, J. Kasahara]

m_dot_H2 = m_dot_tot*H2_Percent;
m_dot_O2 = m_dot_tot*O2_Percent;

disp([' '])
disp(['................................................................']);
disp(['Mass Flow'])

disp([' '])
disp(['Fill Height: ', num2str(Fill_Height),' (m)']);
disp(['Fill Volume: ', num2str(Fill_Volume),' (m^3)']);
disp(['m_dot Total: ', num2str(m_dot_tot),' (kg/s)']);
disp(['m_dot H2: ', num2str(m_dot_H2),' (kg/s)']);
disp(['m_dot O2: ', num2str(m_dot_O2),' (kg/s)']);

%% Thrust Calculation (Check)
% Equations taken from:
% - "Analytical Models for the Thrust of a Rotating Detonation Engine" [J. Shepherd, J. Kasahara]

V_w = CJ_Point(1,1)/(((Minimum_Channel_OD/1000)+(Minimum_Channel_ID/1000))/2);
f = V_w/(2*pi);
T = 1/f;
P_a = 101325;

q_hc = -241800*m_dot_H2*Detonation_Structure(1,23)/0.00202; %heat released based on H2 heat of combustion and the amount of H2 being used
% q_hc = 500000; %q_hc test value for thrust dependancy check

T_t = 2000; %N

% Term_1a = q_hc/(cp1*T1);
% Term_2a = (CJ_Point(1,2)/P1).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
% Term_3a = (P1/CJ_Point(1,2)).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
% Term_4a = (CJ_Point(1,3)/T1);
% T_ue = m_dot_tot*sqrt(2*cp1*T1)*((1 + Term_1a - Term_2a * Term_3a * Term_4a).^0.5);
%assume CJ pressure is exhaust pressure due to short lenght of engine and no nozzle design

Term_1b = q_hc/(cp1*T1)
Term_2b = (P_a/P1).^((CJ_Point(1,14)-1)/CJ_Point(1,14))
Term_3b = (P1/CJ_Point(1,2)).^((CJ_Point(1,14)-1)/CJ_Point(1,14))
Term_4b = (CJ_Point(1,3)/T1)
T_e = m_dot_tot*sqrt(2*cp1*T1)*((1 + Term_1b - Term_2b * Term_3b * Term_4b).^0.5);

% Tsp_ue = T_ue/m_dot_tot;
Tsp_e = T_e/m_dot_tot;

disp([' '])
disp(['................................................................']);
disp(['Thrust'])

disp(['Angular Velocity: ', num2str(V_w),' (rad/s)']);
disp(['Wave Frequency: ', num2str(f),' (Hz)']);
disp([' '])
disp(['Target Thrus: ', num2str(T_t),' (N)']);
% disp(['Actual Thrust (Under Expanded): ', num2str(T_ue), ' (N)']);
disp(['Actual Thrust (Expanded): ', num2str(T_e), ' (N)']);
% disp(['Specific Thrust (Under Expanded): ', num2str(Tsp_ue), ' (N/kg/s)']);
disp(['Specific Thrust (Expanded): ', num2str(Tsp_e), ' (N/kg/s)']);

%% Isp Calculation
%Equations taken from:
% - "SDToolbox: Numerical Tools for Shock and Detonation Wave Modeling" 
% [S. Kao, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson, J. E. Shepherd]

% Isp_ue = T_ue/(m_dot_tot*9.81);

Isp_e = T_e/(m_dot_tot*9.81);

disp([' '])
disp(['................................................................']);
disp(['Specific Impulse'])

disp([' '])
% disp(['Isp (Under Expanded): ', num2str(Isp_ue),' (s)']);
disp(['Isp (Expanded): ', num2str(Isp_e),' (s)']);

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
