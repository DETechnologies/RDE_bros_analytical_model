% The new analytical model!
% DETechnologies!
% Logan and Shak - 2023/2024

close all force 
clear 
clc

mechFiles={'Burke2012.yaml';'h2o2.yaml';'Hong2011.yaml'};
mech = mechFiles{1}; %1 means it uses the first one (burke) (first row (row matrix))

Pressure_range=[180e+3,8*101.325e+3,1e+3]; % low,high,step size [Pa]
Temp_range=[300,300,1]; % low,high,step size [K]
eqv_ratio_range=[1,1,0.05];  % low,high,step size 
CellSizeCorrelationIndex=4;
GeometryCorrelationIndex=1;


CellSizeCorrelations={'Gavrikov','Westbrook','Ng','SeanCB'}; %[1-4]
GeometryCorrelations={'Ng','Bykovskii'}; %[0-1]
n=0;
Outputs=array2table(zeros(0,15),'VariableNames',{'I/P Pressure (Pa)','I/P Temperature (K)','Eqv Ratio','I/P Density (kg/m^3)','Speed of Sound in Propellant (m/s)',...
                                        'CJ Speed (m/s)','VN Pressure (Pa)','CJ Temperature (K)','CJ Pressure (Pa)','Chosen Cell Size Correlation',...
                                        'Cell Size value','WaveNumber','Thrust O/P','ISP','mDot'});
for P1=Pressure_range(1,1):Pressure_range(1,3):Pressure_range(1,2)
    for T1=Temp_range(1,1):Temp_range(1,3):Temp_range(1,2)
        for eq=eqv_ratio_range(1,1):eqv_ratio_range(1,3):eqv_ratio_range(1,2)
            n=n+1;
            [gas1,VN,CJ,ZND,CellSizePredictions,Misc,GeometryPredictor] = NewAnalyticalModel(P1,T1,eq,mech,0,CellSizeCorrelationIndex,GeometryCorrelationIndex,0);
            
            %make a table of the things we actually care about here.
            current=[P1,T1,eq,density(gas1),soundspeed_fr(gas1),CJ(1),VN(2),CJ(3),CJ(2),{sprintf('%s using %s',CellSizeCorrelations{CellSizeCorrelationIndex},GeometryCorrelations{GeometryCorrelationIndex+1})},...
                     GeometryPredictor{CellSizeCorrelationIndex*2-GeometryCorrelationIndex,'CellSize'},Misc{1,'Wave_Number_Sean'},Misc{1,'Thrust'},Misc{1,"ISP"},GeometryPredictor{CellSizeCorrelationIndex*2-GeometryCorrelationIndex,'Mass Flow Rate kg/s'}];

            Outputs=[Outputs;current];
            if ~ mod(n,10) %save o/p every n loops so its fast, but that we dont lose too much data when errors.
                save('AnalyticalModel_calculatorOutput.mat','Outputs') 
            end
        end
    end 
end 