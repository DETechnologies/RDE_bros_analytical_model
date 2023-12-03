%% Function Definition
function Results = ZND_Structure_Shak(P1, T1, q, mech, gas1)
%% SET SHOCK SPEED
%[cj_speed, curve, ~, dnew, plot_data] = CJspeed(P1, T1, q, mech);
[cj_speed] = CJspeed(P1, T1, q, mech);
% CJspeed_plot(2,plot_data,curve,dnew)

%%
% FIND EQUILIBRIUM POST-SHOCK STATE FOR GIVEN SPEED
[gas] = PostShock_eq(cj_speed, P1, T1, q, mech);
u_cj = cj_speed*density(gas1)/density(gas);

% FIND FROZEN PRE-SHOCK STATE FOR GIVEN SPEED
[gas] = PostShock_fr(cj_speed, P1, T1, q, mech);

%% SOLVE ZND DETONATION ODES
[out] = zndsolve(gas,gas1,cj_speed,'advanced_output',true,'t_end',2e-3);

%% Find CV parameters including effective activation energy
set(gas,'Temperature',T1,'Pressure',P1,'X',q);
gas = PostShock_fr(cj_speed, P1, T1, q, mech);
Ts = temperature(gas); 
Ps = pressure(gas);
Ta = Ts*(1.02);
set(gas, 'T', Ta, 'P', Ps, 'X', q);
[CVout1] = cvsolve(gas);
Tb = Ts*(0.98);
set(gas, 'T', Tb, 'P', Ps, 'X', q);
[CVout2] = cvsolve(gas);

%% Approximate effective activation energy for CV explosion
taua = CVout1.ind_time;
taub = CVout2.ind_time;
if(taua==0 || taub==0)
    theta_effective_CV = 0;
else
    theta_effective_CV = 1/Ts*((log(taua)-log(taub))/((1/Ta)-(1/Tb)));
end

%%
%  Find Gavrikov induction length based on 50% limiting species consumption,
%       fuel for lean mixtures, oxygen for rich mixtures
%  Westbrook time based on 50% temperature rise
limit_species = 'H2';
i_limit = speciesIndex(gas,limit_species);
set(gas,'Temperature',Ts,'Pressure',Ps,'X',q);
X = moleFractions(gas);
X_initial = X(i_limit);
equilibrate(gas,'UV');
X = moleFractions(gas);
X_final = X(i_limit);
T_final = temperature(gas);
X_gav = 0.5*(X_initial - X_final) + X_final;
T_west = 0.5*(T_final - Ts) + Ts;
b = length(CVout1.speciesX(:,i_limit));
for i = 1:b
    if (CVout1.speciesX(i:i,i_limit) > X_gav)
        t_gav = CVout1.time(i);
    end
end
x_gav = t_gav*out.U(1);
for i = 1:b
    if (CVout1.T(i) < T_west)
        t_west = CVout1.time(i);
    end
end
x_west = t_west*out.U(1);

max_thermicity_width_ZND = u_cj/out.max_thermicity_ZND;   %Ng et al definition
chi_ng = theta_effective_CV*out.ind_len_ZND/max_thermicity_width_ZND;
cell_gav = gavrikov(x_gav,theta_effective_CV, Ts, T1);
cell_ng = ng(out.ind_len_ZND, chi_ng);

b = length(out.T);

%%  Generate Results Array

Results = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

Results([1]) = Ts;
Results([2]) = Ps;
Results([3]) = Ta; 
Results([4]) = Tb; 
Results([5]) = taua;
Results([6]) = taub;
Results([7]) = theta_effective_CV;
%Results([8]) = limit_species;
%Results([9]) = i_limit; 
Results([10]) = X_initial;
Results([11]) = X;
Results([12]) = b;
Results([13]) = X_final;
Results([14]) = T_final;
Results([15]) = X_gav;
Results([16]) = T_west;
Results([17]) = x_gav;
Results([18]) = x_west;
Results([19]) = max_thermicity_width_ZND;
Results([20]) = chi_ng;
Results([21]) = cell_ng;
Results([22]) = cell_gav;
Results([23]) = out.exo_time_ZND;

%% Display Results

disp([' ']);
disp(['................................................................']);
disp(['ZND computation results ']);

disp([' ']);
disp(['T (K): ??? = ' num2str(out.T(1),5)... 
    newline '       final = ' num2str(out.T(b),5)... 
    newline '       max = ' num2str(max(out.T(:)),5)]);
disp(['P (Pa): ??? = ' num2str(out.P(1),3)...
    newline '        final = ' num2str(out.P(b),3)...
    newline '        max = ' num2str(max(out.P(:)),3)]);
disp(['M: ??? = ' num2str(out.M(1),3)...
    newline '   final = ' num2str(out.M(b),3)...
    newline '   max = ' num2str(max(out.M(:)),3)]);
disp(['u (m/s): ??? = ' num2str(out.U(1),5)...
    newline '         final = ' num2str(out.U(b),5)...
    newline '         cj = ' num2str(u_cj,5)]);

disp([' ']);
disp(['Reaction zone induction length = ',num2str(out.ind_len_ZND,'%8.3e'),' m']);
disp(['Reaction zone induction time = ',num2str(out.ind_time_ZND,'%8.3e'),' s']);
disp(['Reaction zone thermicity  pulse width (exothermic length) = ',num2str(out.exo_len_ZND,'%8.3e'),' m']);
disp(['Reaction zone thermicity  pulse time (exothermic time) = ',num2str(out.exo_time_ZND, '%8.3e'),' s']);
disp(['Reaction zone width (u_cj/sigmadot_max) = ',num2str(max_thermicity_width_ZND,'%8.3e'),' m']);

disp([' ']);
disp(['CV computation results; ']);
%disp(['Time to dT/dt_max = ',num2str(CVout1.ind_time, '%8.3e'),' s']);
%disp(['Distance to dT/dt_max = ',num2str(CVout1.ind_time*out.U(1), '%8.3e'),' m']);
disp(['Reduced activation energy) = ',num2str(theta_effective_CV,'%8.3e')]);
disp(['Time to 50% consumption = ',num2str(t_gav, '%8.3e'),' s']);
disp(['Distance to 50% consumption = ',num2str(x_gav, '%8.3e'),' m']);
disp(['Time to 50% temperature rise = ',num2str(t_west, '%8.3e'),' s']);
disp(['Distance to 50% temperature rise = ',num2str(x_west, '%8.3e'),' m']);

disp([' ']);
disp(['Cell size predictions ']);
disp(['Westbrook correlation ',num2str(29*x_west, '%8.3e'),' m']);
disp(['Cell width estimate (A=29) ',num2str(29*out.ind_len_ZND, '%8.3e'),' m']);
%disp(['Cell size (lambda)',mun2str(lambda), '%8.3e']);
disp(['Gavrikov correlation ',num2str(cell_gav, '%8.3e'),' m']);
disp(['Ng et al Chi Parameter (not cell size!!!) ',num2str(chi_ng, '%8.3e'),' m']);
disp(['Ng et al correlation ',num2str(cell_ng, '%8.3e'),' m']);
% 
% znd_plot(out,'maxx',0.1,'major_species',{'H2', 'O2', 'H2O'},...
%     'minor_species',{'H', 'O', 'OH', 'H2O2', 'HO2'});

%%
function lambda = gavrikov(delta,theta, Tvn, T0)
% Correlation function for detonation cell width 
% proposed by Gavrikov et al COMBUSTION AND FLAME 120:19-33 (2000)
% based on using a reaction zone length based on time to 50% limiting
%       reactant consumption in constant volume explosion approximation using vn
%       postshock velocity to convert time to distance.   Tested against a range
%       of fuel-oxidizer diluent mixtures
%
% Inputs:
%   delta = reaction zone length based on time to 50% consumption of limiting
%   reactant from CV computation and delta = time * w_VN
%   theta = Ea/RT_VN,  effective reduced activation energy based on CV computation
%   Tvn = von Neumann (postshock temperature behind CJ shock wave)
%   T0 = initial temperature
%
% Constants
a = -0.007843787493;
b = 0.1777662961;
c = 0.02371845901;
d = 1.477047968;
e = 0.1545112957;
f = 0.01547021569;
g = -1.446582357;
h = 8.730494354;
i = 4.599907939;
j = 7.443410379;
k = 0.4058325462;
m = 1.453392165;
%  define nondimensional parameters
X = theta;
Y = Tvn/T0;
z = Y*(a*Y-b) + X*(c*X-d + (e-f*Y)*Y) + g*log(Y)+ h*log(X) + Y*(i/X - k*Y/X^m) - j;
lambda = delta*10^(z);
end

%%
function lambda = ng(delta,chi)
% correlation function for detonation cell size from
%    Ng, Hoi Dick, Yiguang Ju, and John H. S. Lee. 2007. Assessment of
%    Detonation Hazards in High-Pressure Hydrogen Storage from Chemical
%   Sensitivity Analysis. INTERNATIONAL JOURNAL OF HYDROGEN ENERGY 32 (1): 93-99.
% Tested only against low pressure H2-air data

% Inputs:
%   delta = reaction zone length based on peak thermicity in ZND simulation
%   chi = theta*Delta_i/Delta_r where 
%       theta = reduced effective activation energy from CV computation
%       Delta_i = distance to peak thermicity from ZND computation
%       Delta_r = w_vN/\sigmadot_max from ZND computation
%   See Ng et al.  Combustion Theory and Modeling 2005 for a discussion of the chi parameter.  
%
% Constants
A0 = 30.465860763763;
a1 = 89.55438805808153;
a2 = -130.792822369483;
a3 = 42.02450507117405;
b1 = -0.02929128383850;
b2 = 1.0263250730647101E-5;
b3 = -1.031921244571857E-9;
lambda = delta*(A0 + ((a3/chi + a2)/chi + a1)/chi + ((b3*chi + b2)*chi + b1)*chi);
end

end
