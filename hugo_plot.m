%% HSKP
close all 
clear 
clc 

%% read data from excel
% E = readtable("C:\Users\Logan\MATLAB Drive\RDE Bros\RDE initial params.xlsx");
% E = readtable("C:\Users\Shakib\MATLAB Drive\RDE Bros\RDE initial params.xlsx"); %Shak Laptop
E = readtable("E:\MatLab Drive\RDE Bros\RDE initial params.xlsx"); %Shak Desktop

%% define variables 
% Define the range of rho2 values for which we want to plot the curve
% rho1=E{4,2} %[kg/m^3] equivalent input density of propellant
rho1=4.8149; %from anal model

% v=E{5,2}; %Specific Heat ratio
v = 1.3108; %from anal model 

% p1=E{1,5}*1000; %[kPa] input pressure
p1 = 1000e+3; %from anal model

% q=E{6,2}; %energy released in combustion
q = 1.372311620796457e+07; %from anal model

rho2 = linspace(0.05, 100, 1000);
% rho2 = 8.8186; %from anal model

%% Define shak's updated function
A=(v/(v-1));
ps=(1./rho2)+(1/rho1);
C=(0/A)+(p1/rho1)-((p1/(2*A)).*ps);
D=(1./rho2)-((1/(2*A)).*ps);
P2=C./D;

c_c=(q/A)+(p1/rho1)-((p1/(2*A)).*ps);
P2_c=c_c./D;

%% Define Logan's alternate formulation for hugoniot plot
% hugo formula works perfectly
% % a=A;
% % b=(1/rho1 +1./rho2);
% % P2_lp=(q+(p1*a/rho1)-(b*p1/2))./((a./rho2)-(b/2));

%%cj fuckers


%% other graphical things:
Pmin=p1+(v-1)*q*rho1;
v_asymp=rho1*((v+1)/(v-1));
h_asymp=-p1*((v-1)/(v+1)); %no sense graphing its negative anyways
%% define CJ points
E=(2*v*p1)/(q*rho1*((v^2)-1));
CJ(1,1)=p1+(v-1)*q*rho1*(1+sqrt(1+E)); %y coord upper
CJ(1,2)=inv(rho1)+((v-1)/v)*(q/p1)*(1-sqrt(1+E)); %X-coord upper

CJ(2,1)=p1+(v-1)*q*rho1*(1-sqrt(1+E)); %y coord lwr
CJ(2,2)=inv(rho1)+((v-1)/v)*(q/p1)*(1+sqrt(1+E)); %X-coord lwr
%% Define Rayleigh Lines
rayleigh_slope=(CJ(1,1)-p1)/(CJ(1,2)-1/rho1);
rayleigh_y_int=p1-rayleigh_slope*(1/rho1);

%% Plot
% Shock Hugoniot curve
figure('Name','Hugoniot Curves');
plot(1./rho2, P2,'color','black'); %this is the shock hugoniot
hold on
% combustion Hugoniot curve
plot(1./rho2, P2_c,'color',[128/255 128/255 128/255]); %this is the combustion hugoniot
plot(1/rho1,p1,'r*'); %input point

ip_callout=sprintf('%.2e', p1)+" Pa,"+sprintf('%.2e', 1/rho1)+"m^3/kg";
% text(1/rho1+0.01,p1+0.01,ip_callout)

%%CJ points
plot(CJ(1,2),CJ(1,1),'b*') %upper
plot(CJ(2,2),CJ(2,1),'black*') %lower

%%p_inf_min as an asymp
x=linspace(0,1.5,2);
% plot(x,Pmin*ones(size(x)),':black');

%%plot rayleigh line
x=linspace(0,1/rho1,2);
y=rayleigh_slope*x+rayleigh_y_int;
plot(x,y,':','Color',[128/255 128/255 128/255]);

%%plot intersection of rayleigh and shock hugo lines.
warning('off')
[xint,yint]=intersections(x,y,(1./rho2),P2);
plot(xint(2,1),yint(2,1),'blackx') 
callout="P_{VN}: "+sprintf('%.2e', yint(2, 1))+" Pa";
% text(xint(2,1)+0.02,yint(2,1),callout)

detonationCallout="P_{det}: "+sprintf('%.2e', CJ(1, 1))+" Pa";
% text(CJ(1,2)+0.01,CJ(1,1),detonationCallout)

%make some y-space gear for the vertical asymptote
y=linspace(0,1.2*10^7,2);
% plot((1/v_asymp)*ones(size(y)),y,'--black');

ylim([0,2.5e+7])
xlim([0.027,0.45])
xlabel('Specific Volume (1/rho) [m^3/kg]');
ylabel('Pressure (P) [Pa]');
title('Hugoniot Curve for Constant Heat Release');
legend('Shock Hugoniot','Combustion Hugoniot','Input','Upper CJ','LowerCJ','Rayleigh Line','VonNeumann Spike')

%% Plot one to print in our paper 
% Shock Hugoniot curve
figure('Name','Hugoniot Curve - Markup with Zones');
% plot(1./rho2, P2,'color','black'); %this is the shock hugoniot
hold on
% % combustion Hugoniot curve
plot(1./rho2, P2_c,'color',[128/255 128/255 128/255]); %this is the combustion hugoniot
% plot(1/rho1,p1,'x',Color='red'); %input point

% ip_callout=sprintf('%.2e', p1)+" Pa,"+sprintf('%.2e', 1/rho1)+"m^3/kg";
% % text(1/rho1+0.01,p1+0.01,ip_callout)

%%CJ points
% plot(CJ(1,2),CJ(1,1),'o',Color='black') %upper
% plot(CJ(2,2),CJ(2,1),'o',Color='black') %lower

% %%p_inf_min as an asymp
% x=linspace(0,1.5,2);
% plot(x,Pmin*ones(size(x)),':black');

% %%plot rayleigh line
% x=linspace(0,1/rho1,2);
% y=rayleigh_slope*x+rayleigh_y_int;
% plot(x,y,'--','Color',[0 0.4470 0.7410]);

%%plot intersection of rayleigh and shock hugo lines.
% warning('off')
% [xint,yint]=intersections(x,y,(1./rho2),P2);
% plot(xint(2,1),yint(2,1),'rx') 
% callout="P_{VN}: "+sprintf('%.2e', yint(2, 1))+" Pa";
% text(xint(2,1)+0.02,yint(2,1),callout)

% detonationCallout="P_{det}: "+sprintf('%.2e', CJ(1, 1))+" Pa";
% text(CJ(1,2)+0.01,CJ(1,1),detonationCallout)

%make some y-space gear for the vertical asymptote
% y=linspace(0,1.2*10^7,2);
% plot((1/v_asymp)*ones(size(y)),y,'--black');

ylim([0.2,2.5e+7])
xlim([0.06,0.45])
xlabel('Specific Volume (1/rho) [m^3/kg]');
ylabel('Pressure (P) [Pa]');
title('Hugoniot Curve for Constant Heat Release');
legend('Combustion Hugoniot','Input','Upper CJ','LowerCJ')
