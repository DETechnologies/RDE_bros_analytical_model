close all force
clear 
clc

%  Matlab DOE
load("Output_data_forDOE.mat")

dbb=bbdesign(3);
ShowAll=0
if ShowAll
    figure("Name","BBD DOE Design")
    plot3(dbb(:,1),dbb(:,2),dbb(:,3),'ro','MarkerFaceColor','b')
    
    X=[1 -1 -1 -1 1 -1 -1 -1 1 1 -1 -1;
        1 1 1 -1 1 1 1 -1 1 1 -1 -1];
    Y=[-1 -1 1 -1 -1 -1 1 -1 1 -1 1 -1;
        1 -1 1 1 1 -1 1 1 1 -1 1 -1];
    
    Z=[1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1;
        1 1 1 1 -1 -1 -1 -1 1 1 1 1];
    
    line(X,Y,Z,'Color','b')
end
%% now do real things

OT=Output';

% figure("Name","Interaction Plot")
y=randn(1000,1);
% g=ceil(3*rand(9,24));
% interactionplot(y,g)
% interactionplot(Output(:,1),g)

% figure("Name","Multivariable chart")
% multivarichart(Output(:,1:2),Output(:,4:5))

% anova1(Output(:,:))
anova1(OT(:,:))
