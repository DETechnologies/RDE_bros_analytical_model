close all force
clear
clc

P1 = 10000; % [Pa]
T1 = 293;% [K]
eq=1.0;
mech = 'Burke2012.yaml';

gas1 = Solution(mech);
eq=InitialState(T1,P1,eq,gas1); %%this calculates the mol ratio of hydrogen to oxygen
FAR = sprintf('H2:%d O2:%d',eq(1,1),eq(2,1));
H2_Percent = eq(1,1);
O2_Percent = eq(2,1);

[cj_speed] = CJspeed(P1, T1, FAR, mech);

