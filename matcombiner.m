pt1=load('Output_data_dec10_new.mat');
pt2=load('Output_data_dec10_new_pt2.mat');
% pt3=load('Output_data_dec10.mat');

Output=[pt1.Output;pt2.Output];
Output_dataNames=pt1.Output_dataNames;
save('Output_data_dec11_new_all.mat','Output',"Output_dataNames")