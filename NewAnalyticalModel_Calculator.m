% The new analytical model
% Logan and Shak - 2023 

%% init params
P1 = 180e+3; % [Pa]
T1 = 300;% [K]
eq = 1.0;
% mech = 'Burke2012.yaml';
mech = 'h2o2.yaml';
% mech = 'Hong2011.yaml';
% mech = 'sandiego20161214_H2only.yaml'; % this is fucked

%% Call function
OP=NewAnalyticalModel(P1,T1,eq,mech);