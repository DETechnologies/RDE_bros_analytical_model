%% Set Input Parameters
clear;
clc;
close all;
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

Output_dataNames=["I/P Pressure [pa]","I/P Temperature [K]","Eqv Ratio","I/P Density"," I/P Gamma","Speed of Sound (in that gas mix) [m/s]"];
sz=0;%initializes the size fcn, so it doesnt error out

%% Data Aq. Loop

for p=Pressure_range(1,1):Pressure_range(1,3):Pressure_range(1,2)
    for t=Temp_range(1,1):Temp_range(1,3):Temp_range(1,2)
        for e=eqv_ratio_range(1,1):eqv_ratio_range(1,3):eqv_ratio_range(1,2)
            InitialState(t,p,e,gas_i) %returns density, 
            gamma=cp_mass(gas_i)/cv_mass(gas_i);
            spd=soundspeed_fr(gas_i);

            Output(sz(1,1)+1,:)=[p,t,e,density(gas_i),gamma,spd];
            sz=size(Output);
        end 
    end
end

%% saves the file as .mat file for later accessing

save('Output_data.mat','Output',"Output_dataNames")