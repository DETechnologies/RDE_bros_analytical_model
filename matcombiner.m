pt1=load('noDelete_mat\OutputVPressure_v3.mat');
pt2=load('noDelete_mat\OutputVPressure_v3_2.mat');
% pt3=load('noDelete_mat\OutputVEqvR_YesSean3.mat');

Output=[pt1.Output;pt2.Output];
Output_dataNames=pt1.Output_dataNames;
save('noDelete_mat\OutputVPressure_v3_combined.mat','Output',"Output_dataNames")