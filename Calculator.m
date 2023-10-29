%% Set Input Parameters
clear
clc
close all
disp('Analytical Model Calculator')

Pressure_range=[1.013e+6,2.533e+6,0.76e+6]; % low,high,step % in pascals
Temp_range=[273.15,473.15,100]; % low,high,step % in Kelvin
eqv_ratio_range=[0.8,1.5,0.1];  % low,high,step %no units
%%% not sure if this matters but make sure the step size is a whole
%%% multiple of the range.

mech = 'h2o2.yaml'; %%yaml files come from here: C:\Program Files\Cantera\data
gas_i = Solution(mech);

%% How data is stored

%%% The Output.mat file is organized as follows
%%% I/P Pressure | I/P Temperature | Eqv. Ratio | I/P density | I/P gamma (heat capacity
%%% ratio) | Speed |
%%% the rows under these columns contain all the data.

Output_dataNames=["I/P Pressure [Pa]","I/P Temperature [K]","Eqv Ratio","I/P Density [kg/m^3]","Speed of Sound (in that gas mix) [m/s]",...
                    "CJ Speed [m/s]","VN Pressure [Pa]","CJ Temperature [K]","CJ Pressure [Pa]","CJ Density [kg/m^3]","Reaction zone induction length [m]",...
                    "Reaction zone induction time [s]","Reaction zone thermicity  pulse width (exothermic length) [m?]","Reaction zone thermicity  pulse time (exothermic time) [s]","Reaction zone width (u_cj/sigmadot_max)",...
                    "Westbrook Cell Size [m]","Cell width estimate (A=29) [m]","Gavrikov Cell Size [m]","Ng et al correlation cell size[m]",...
                    "Final specific impulse (fr) (Isp) [s]"];
sz=0;%initializes the size fcn, so it doesnt error out

%% Data Aq. Loop

for p=Pressure_range(1,1):Pressure_range(1,3):Pressure_range(1,2)
    for t=Temp_range(1,1):Temp_range(1,3):Temp_range(1,2)
        for e=eqv_ratio_range(1,1):eqv_ratio_range(1,3):eqv_ratio_range(1,2)
            tic %start timer
            
            %this section gets initial parameters
            FAR=InitialState(t,p,e,gas_i); %returns adjusted FuelAirRatio (H2:O2) based on eqv ratio 
            spd_sound_gas=soundspeed_fr(gas_i);

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

            %this section gets VN things
            ZNDResults=ZND_Structure(p,t,FAR,mech,gas_i);
            ZND_out=zndsolve(VN_gas,gas_i,CJ_spd,'advanced_output',true,'t_end',2e-3);
%             u_cj=CJ_spd*density(gas_i)/CJ_density;
%             max_thermicity_width_ZND=u_cj/ZND_out.max_thermicity_ZND; %Ng et al definition
            x_west=ZNDResults(18);
            max_thermicity_width_ZND=ZNDResults(19);
            cell_gav=ZNDResults(22);
            cell_ng=ZNDResults(21);

            %this is rocket impulse things
            RResults=Rocket_Impulse(t,p,FAR,mech);
            frozen_Isp=RResults(9);
            
            %this is where all the data gets stored each iteration
            Output(sz(1,1)+1,:)=[p,t,e,density(gas_i),spd_sound_gas,CJ_spd,VN_pressure,...
                                CJ_temp,CJ_pressure,CJ_density,ZND_out.ind_len_ZND,ZND_out.ind_time_ZND,...
                                ZND_out.exo_len_ZND,ZND_out.exo_time_ZND,max_thermicity_width_ZND,...
                                29*x_west,29*ZND_out.ind_len_ZND,cell_gav,cell_ng,frozen_Isp];
            sz=size(Output);
            fprintf('looptime: %d\n',toc)
        end 
    end
end

%% saves the file as .mat file for later accessing

save('Output_data.mat','Output',"Output_dataNames")