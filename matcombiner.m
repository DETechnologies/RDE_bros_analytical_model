wrm=load('Output_data_combined_warm.mat');
cold=load('Output_data_combined_cold.mat');

Output=[cold.Output;wrm.Output];
Output_dataNames=wrm.Output_dataNames;
save('Output_data_combined_all.mat','Output',"Output_dataNames")