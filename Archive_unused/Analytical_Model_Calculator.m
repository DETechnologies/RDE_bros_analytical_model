%% Set Input Parameters
clear
% clc
close all all force
disp('Analytical Model Calculator')

Pressure_range=[500e+3,500e+3,1e+3]; % low,high,step size % in pascals (range 5 atm to 25 atm (500psi))
Temp_range=[300,300,1]; % low,high,step size % in Kelvin
eqv_ratio_range=[1,1,0.05];  % low,high,step size %no units

P_a = 101325;

mech = 'h2o2.yaml'; %%yaml files come from here: C:\Program Files\Cantera\data
gas_i = Solution(mech);
g=9.81;%m/s^2
q_h = (241820/18.01528)*1000;

%% time collection estimate; 
p_steps=(Pressure_range(1,2)-Pressure_range(1,1))/Pressure_range(1,3)+1;
t_steps=(Temp_range(1,2)-Temp_range(1,1))/Temp_range(1,3)+1;
eqv_steps=(eqv_ratio_range(1,2)-eqv_ratio_range(1,1))/eqv_ratio_range(1,3)+1;

total_loops=p_steps*t_steps*eqv_steps;
time_est=(total_loops*10)/60; %takes ~10s per loop

fprintf("%d Number of Loops, Estimated ~ %d mins to complete",total_loops,time_est)
fprintf('\n')

%% How data is stored

%%% The Output.mat file is organized as follows
%%% I/P Pressure | I/P Temperature | Eqv. Ratio | I/P density | I/P gamma (heat capacity
%%% ratio) | Speed |
%%% the rows under these columns contain all the data.

Output_dataNames=["I/P Pressure [Pa]","I/P Temperature [K]","Eqv Ratio","I/P Density [kg/m^3]","Speed of Sound (in that gas mix) [m/s]",...
                    "CJ Speed [m/s]","VN Pressure [Pa]","CJ Temperature [K]","CJ Pressure [Pa]","CJ Density [kg/m^3]","Reaction zone induction length [m]",...
                    "Reaction zone induction time [s]","Reaction zone thermicity  pulse width (exothermic length) [m?]","Reaction zone thermicity  pulse time (exothermic time) [s]","Reaction zone width (u_cj/sigmadot_max)",...
                    "Westbrook Cell Size [m]","Cell width estimate (A=29) [m]","Gavrikov Cell Size [m]","Ng et al correlation cell size [m]",...
                    "Minimum Channel OD [mm]", "Minimum Channel Width [mm]", "Minimum Chamber Length [mm]"...
                    "big_red_minimumFillHeight [mm]", "big_red_minD [mm]", "big_red_min_delta [mm]", "big_red_minLength [mm]", "Fill_Height [m]",...
                    "Fill_Volume [m^3]", "T_e [N]",...
                    "Tsp_e [N/kg/s]", ...
                    "ISP_e [s]","Sean CB Cell Size","sean wave number","sean fill time","Shepherd M_{dot}"];
sz=0;%initializes the size fcn, so it doesnt error out

%% Data Aq. Loop

n=0;
for p=Pressure_range(1,1):Pressure_range(1,3):Pressure_range(1,2)
    for t=Temp_range(1,1):Temp_range(1,3):Temp_range(1,2)
        for e=eqv_ratio_range(1,1):eqv_ratio_range(1,3):eqv_ratio_range(1,2)
            tic %start timer
            n=n+1;

            %this section gets initial parameters
            eq=InitialState(t,p,e,gas_i); %returns adjusted FuelAirRatio (H2:O2) based on eqv ratio 
            FAR = sprintf('H2:%d O2:%d',eq(1,1),eq(2,1));
            spd_sound_gas=soundspeed_fr(gas_i);
%             R1 = density(gas_i);
            gas_downstr=PostShock_eq(spd_sound_gas,p, t, FAR, mech);
            R2=density(gas_downstr);
            cp1 = cp_mass(gas_i);
            
            %this section is an alternative to using the vN_state.m fcn
            %(less variables)
            CJ_spd=CJspeed(p,t,FAR,mech);
            VN_gas=PostShock_fr(CJ_spd,p,t,FAR,mech);
            VN_pressure=pressure(VN_gas);
            
            %this section gets CJ things using alternative to CJ_state fcn
            CJ_gas=PostShock_eq(CJ_spd,p,t,FAR,mech);           
            CJ_temp=temperature(CJ_gas);
            CJ_pressure=pressure(CJ_gas);
            CJ_density=density(CJ_gas);

            %this section gets ZND things
            ZNDResults=ZND_Structure(p,t,FAR,mech,gas_i);
            ZND_out=zndsolve(VN_gas,gas_i,CJ_spd,'advanced_output',true,'t_end',2e-3);
            %             u_cj=CJ_spd*density(gas_i)/CJ_density;
            %             max_thermicity_width_ZND=u_cj/ZND_out.max_thermicity_ZND; %Ng et al definition
            x_west=ZNDResults(18);
            max_thermicity_width_ZND=ZNDResults(19);
            cell_gav=ZNDResults(22);
            cell_ng=ZNDResults(21);
            
            %sean CB stuff
            cell_sean = (1.6e-3*101325)/p;

            %geometry defns
            Minimum_Channel_OD=40*cell_sean;
            Minimum_Channel_Width=2.4*cell_sean;
            Minimum_Channel_ID = (Minimum_Channel_OD - 2*Minimum_Channel_Width);
            Minimum_Chamber_Length=24*cell_sean;

            %geometry big red 
            big_red_minimumFillHeight=(12+0)*cell_gav;
            big_red_minD=28*cell_gav;
            big_red_min_delta=0.2*big_red_minimumFillHeight;
            big_red_minLength=2*big_red_minimumFillHeight;
            
            %m_dots
            Fill_Height = (12-0)*cell_sean; %This is the max critical fill height case
            Fill_Volume = 0.25*pi*(((Minimum_Channel_OD).^2)-((Minimum_Channel_ID).^2))*Fill_Height;
%             m_dot_tot_exit = Fill_Height*(Minimum_Channel_Width/1000)*CJ_density*CJ_spd; %density of combustion products, and cj speed
%             C=pi()*(Minimum_Channel_OD/1000);
%             m_dot_intoengine=(Fill_Volume*R2*CJ_spd/C);
            
            % not proud of this (ignore)
            CJ_Point = CJ_State(p, t, FAR, mech, gas_i,1);
            P1=p;
            T1=t;

            m_dot_P_history = Fill_Height*(Minimum_Channel_Width)*density(gas_i)*CJ_Point(1,1); %density of combustion products, and cj speed
            % [J. Shepherd, J. Kasahara]

            % thrust
            Term_1b = q_h/(cp1*T1);
            Term_2b = (P_a/P1).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
            Term_3b = (P1/CJ_Point(1,2)).^((CJ_Point(1,14)-1)/CJ_Point(1,14));
            Term_4b = (CJ_Point(1,3)/T1);
            T_e = m_dot_P_history*(sqrt(2*cp1*T1))*(sqrt(1 + Term_1b - Term_2b * Term_3b * Term_4b));
            Tsp_e = T_e/m_dot_P_history;

            %spec impulse
            Isp_e=T_e/(m_dot_P_history*9.81);
            
            %more sean cb things
            R_sp = 8.314462618;
            Cl=12; % this is the 12+-5 thing
            m_dot_P_history = Fill_Height*(Minimum_Channel_Width/1000)*CJ_density*CJ_spd; %density of combustion products, and cj speed % [J. Shepherd, J. Kasahara]
            Wave_Number_Sean = (m_dot_P_history*R_sp*T1)/(Cl*0.0016*101325*CJ_spd*Fill_Height);
            Mean_Channel_Diam=Minimum_Channel_OD-Minimum_Channel_Width/2;
            Fill_Time_Sean = (pi*(Mean_Channel_Diam))/(CJ_spd*Wave_Number_Sean);


            %this is where all the data gets stored each iteration
            Output(sz(1,1)+1,:)=[p,t,e,density(gas_i),spd_sound_gas,CJ_spd,VN_pressure,...
                                CJ_temp,CJ_pressure,CJ_density,ZND_out.ind_len_ZND,ZND_out.ind_time_ZND,...
                                ZND_out.exo_len_ZND,ZND_out.exo_time_ZND,max_thermicity_width_ZND,...
                                29*x_west,29*ZND_out.ind_len_ZND,cell_gav,cell_ng...
                                Minimum_Channel_OD,Minimum_Channel_Width,Minimum_Chamber_Length...
                                big_red_minimumFillHeight,big_red_minD,big_red_min_delta,big_red_minLength...
                                Fill_Height,Fill_Volume,...
                                T_e,Tsp_e,Isp_e,cell_sean,Wave_Number_Sean,Fill_Time_Sean,m_dot_P_history];

            sz=size(Output);
            fprintf('\n \nlooptime: %d, loop number: %d of %d',toc,n,total_loops)
            save('Output_troubleshooting.mat','Output',"Output_dataNames") %% saves the file as .mat file for later accessing saving here incase it errors out, overwrite each loop
        end 
    end
end
