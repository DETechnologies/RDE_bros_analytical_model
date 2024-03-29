%% Function Definition
function Results = CJ_State(P1, T1, q, mech, gas1,PrintThings)

%% Find CJ speed
[cj_speed] = CJspeed(P1, T1, q, mech);

%% Evaluate downstream gas state
[gas] = PostShock_eq(cj_speed,P1, T1, q, mech);

%% Evaluate downstream properties of gas object 
Results = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];

T2 = temperature(gas);
P2 = pressure(gas);
R2 = density(gas);
V2 = 1/R2;
S2 = entropy_mass(gas);
H2 = enthalpy_mass(gas);
w2 = density(gas1)*cj_speed/R2;
u2 = cj_speed - w2;
x2 = moleFractions(gas);
c2_eq = soundspeed_eq(gas);
c2_fr = soundspeed_fr(gas);
gamma2_fr = c2_fr*c2_fr*R2/P2;
gamma2_eq = c2_eq*c2_eq*R2/P2;
MCJ = cj_speed/c2_eq;

Results([1]) = cj_speed;
Results([2]) = P2;
Results([3]) = T2;
Results([4]) = R2;
Results([5]) = V2;
Results([6]) = S2;
Results([7]) = H2;
Results([8]) = w2;
Results([9]) = u2;
%Results([10]) = x2;
Results([11]) = c2_eq; %(speed of sound of the product mixture)
Results([12]) = c2_fr; %(speed of sound of the product mixture)
Results([13]) = gamma2_fr;
Results([14]) = gamma2_eq;
Results([15]) = MCJ;

%% Print out
if PrintThings
    disp([' '])
    disp(['................................................................']);
    disp( 'CJ Point Properties');
    % 
    % disp([' '])
    disp(['   CJ speed: ',num2str(cj_speed),' (m/s)']);
    disp(['   Pressure: ',num2str(P2),' (Pa)']);
    disp(['   Temperature: ',num2str(T2),' (K)']);
    disp(['   Density: ',num2str(R2),' (kg/m3)']);
    %disp(['   Entropy: ',num2str(S2),' (J/kg-K)']);
    disp(['   Enthalpy: ',num2str(H2),' (J)']);
    disp(['   Mole Fractions: H2:',num2str(x2(1,1)),' & O2:',num2str(x2(4,1)), '(mol/mol)'])
    % disp(['   w2 (wave frame): ',num2str(w2),' (m/s)']);
    % disp(['   u2 (lab frame): ',num2str(u2),' (m/s)']);
    % disp(['   c2 (frozen) (speed of sound of the mixture): ',num2str(c2_fr),' (m/s)']);
    disp(['   c2 (equilibrium) (speed of sound for mixture): ',num2str(c2_eq),' (m/s)']);
    % disp(['   gamma2 (frozen): ',num2str(gamma2_fr)]);
    disp(['   gamma2 (equilibrium): ',num2str(gamma2_eq)]);
end 

end