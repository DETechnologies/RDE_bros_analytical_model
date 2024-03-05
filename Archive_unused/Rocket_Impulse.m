%% Function Definition
function Results = Rocket_Impulse(T1, P1, q, mech)

%% Set Gas States
gas = Solution(mech);
gas1 = Solution(mech);
set(gas, 'T', T1, 'P', P1, 'X', q);

%% Compute equilibrium properties for constant pressure burn
equilibrate(gas,'HP');
qc = moleFractions(gas);

%% Isentropic expansion from stagnation state to low pressure
Results = [0,0,0,0,0,0,0];

Tt = temperature(gas);
St = entropy_mass(gas);
Ht = enthalpy_mass(gas);
Pt = pressure(gas);
Rt = density(gas);
ct = soundspeed_eq(gas);
gammat_eq = ct^2*Rt/Pt;

Results([1]) = Pt;
Results([2]) = Tt;
Results([3]) = Rt;
Results([4]) = Ht;
Results([5]) = St;
Results([6]) = ct;
Results([7]) = gammat_eq;

% %% Display Ideal Stagnation Properties
% disp(' ');
% disp('-------------------------------------------');
% disp( 'Ideal Rocket Impulse');
% 
% disp(' ');
% disp(['Total pressure ',num2str(Pt),' (Pa)']);
% disp(['Total temperature ',num2str(Tt),' (K)']);
% disp(['Total density ',num2str(Rt),' (kg/m3)']);
% disp(['Total enthalpy ',num2str(Ht),' (J/kg-K)']);
% %disp(['Total entropy ',num2str(St),' (J/kg-K)']);
% disp(['Total sound speed (equilibrium) ',num2str(ct),' (m/s)']);
% disp(['gamma2 (equilibrium) ',num2str(gammat_eq),' (m/s)']);

%% ambient pressure for impulse computation
Pa = 0.0;

%%  Isentropic Expansion Computation
pp = Pt;
imax = 120;
for i = 1:1:imax
     pp = pp*0.95;
     sp = [St, pp];
     
    % compute equilibrium isentrope
    setState_SP(gas,sp);
    
    %equilibrate(gas,'SP','maxsteps 10000','maxiter 10000');
    equilibrate(gas,'SP');
    P(i) = pressure(gas);
    R(i) = density(gas);
    T(i) = temperature(gas);
    h_eq(i) = enthalpy_mass(gas);
    u_eq(i) = sqrt(2*(Ht-h_eq(i)));
    Isp_eq(i) = (u_eq(i)+ (P(i) - Pa)/(u_eq(i)*R(i)))/9.81 ; 
    
    % compute frozen impulse (first call is so we get a good guess for entropy 
    set(gas1,'Pressure',pp,'Temperature',T(i),'MoleFractions',qc);
    set(gas1,'Pressure',pp,'Entropy',St,'MoleFractions',qc);
    h_fr(i) = enthalpy_mass(gas1);
    R_fr(i) = density(gas1);
    T_fr(i) = temperature(gas1);
    u_fr(i) = sqrt(2*(Ht-h_fr(i)));
    Isp_fr(i) = (u_fr(i)+ (P(i) - Pa)/(u_fr(i)*R_fr(i)))/9.81 ; 
    
    
end

Results([8]) = Isp_eq(imax);
Results([9]) = Isp_fr(imax);
Results([10]) = P(imax);
Results([11]) = u_eq(imax);
Results([12]) = u_fr(imax);

% %% Display Expanded Properties
% disp(' ');
% disp(['final pressure ',num2str(P(imax)),' (Pa)']);
% disp(['final specific impulse (eq) ',num2str(Isp_eq(imax)),' (s)']);
% disp(['final velocity (eq) ',num2str(u_eq(imax)),' (m/s)']);
% disp(['final specific impulse (fr) ',num2str(Isp_fr(imax)),' (s)']);
% disp(['final velocity (fr) ',num2str(u_fr(imax)),' (m/s)']);