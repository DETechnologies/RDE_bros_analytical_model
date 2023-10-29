clear;
clc;
close all;
disp('3D_Plot_Test')

mech = 'h2o2.yaml';

X = [];
Y = [];
Z = [];
Outputs = [];

H2_s = 2;
O2_s = 1;
FAR_s = (H2_s/O2_s);

% keep* length of iteration intervals the same
for P1 = 500e3:100e3:1000e3
    X = [X, P1];
    
    for T1 = 200:10:300
        Y = [Y, T1];
        
        for eqv = 0.1:0.4:2
            Z = [Z, eqv];

            FAR_a = eqv*FAR_s;
            H2_a = FAR_a;
            O2_a = O2_s;
            q = 'H2:H2_a, O2:O2_a';
            
            gas1 = Solution(mech);
            set(gas1,'Temperature',T1,'Pressure',P1,'MoleFractions',q);
            
            CJ_Point = CJ_State(P1, T1, q, mech, gas1);
            
            Outputs = [Outputs, CJ_Point];
            Outputs = [Outputs, newline];

            Z = [Z, newline];

        end
    
    Y = [Y, newline];

    end

X = [X, newline];

end
