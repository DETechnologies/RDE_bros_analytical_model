%% Function Definitions -- used in the calculator
function FAR = InitialState(T1,P1,e1,gas1)

q = [2;1]; %exact stoichiometric molar ratio is 2/3;1/3 (Hydrogen;oxygen)

h=q(1,1)*e1; %hydrogen as a molar ratio
eq(1,1)=h/(h+1); %hydrogen being converted to percentage
eq(2,1)=1/(h+1); %oxygen mols being converted to percentage

FAR= sprintf('H2:%d O2:%d',eq(1,1),eq(2,1)); %FuelAirRatio as percentage

% disp(eqv_ratio)

%this section is used to prove that the right equivalence ratio is being
% used; we divide by 2 to account for the orignal stoic ratio of the mixture
% r=eq(1,1)/eq(2,1);
% disp(r/2)

set(gas1,'Temperature',T1,'Pressure',P1,'MoleFractions',FAR);