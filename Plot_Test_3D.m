clear;
clc;
close all;
disp('3D_Plot_Test')

mech = 'h2o2.yaml';

X = [];
Y = [];
Z = [];

H2_s = 2;
O2_s = 1;
FAR_s = (H2_s/O2_s);

% keep* length of iteration intervals the same
for i = 500000:100000:1000000
    X = [X, i];
    
    for j = 200:10:300
        Y = [Y, j];
        
        for k = 0.1:0.4:2
            Z = [Z, k];

            FAR_a = k*FAR_s;
            H2_a = FAR_a;
            O2_a = O2_s;
            q = 'H2:H2_a, O2:O2_a';
            
            gasi = Solution(mech);
            set(gasi,'Temperature',j,'Pressure',i,'MoleFractions',q);
            
            CJ_Point = CJ_State(i, j, q, mech, gasi);
            
            append results to matrix

            Z = [Z, newline];

        end
    
    Y = [Y, newline];

    end

X = [X, newline];

end
