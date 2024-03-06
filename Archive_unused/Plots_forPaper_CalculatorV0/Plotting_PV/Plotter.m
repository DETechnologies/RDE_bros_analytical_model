close all 
clear
clc 

makePlots=1

%% PV
PV_raw=readtable("thermo_data.csv",Delimiter=',');

PV=table2array(PV_raw(:,1:2));
PV_labels=table2array(PV_raw(:,3));

p=unique(PV_labels);
a=find(strcmp(PV_labels,p(1,1)));

if makePlots==1
    figure(Name="PV")
end 

masterData=[];
for i = 1:length(p)
    t=find(strcmp(PV_labels,p(i)));
    
    vals=PV(min(t):max(t),:);
    assignin('caller',char(p(i)),cell2table(num2cell(vals),'VariableNames',{'Pressure','Volume',}))
%     assignin('caller',char(p(i)),vals)
    labels(i,1)=p(i);

    if makePlots==1
        semilogy(vals(:,1),vals(:,2))
        grid on
        hold on
    end 
end 

if makePlots==1
%     legend(labels(1:size(labels,1)))
    legend("Brayton (Isobaric)","Fickett-Jacobs (Detonation)","Humphrey (Isochoric)")
    xlabel('Volume [m^3]')
    ylabel('Pressure [bar] ')
    set(gca, 'YTick', [5,10,20,30,40,50])
end 


