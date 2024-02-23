pt1=load('noDelete_mat\OutputVEqvR_YesSean.mat');
pt2=load('noDelete_mat\OutputVEqvR_YesSean2.mat');
pt3=load('noDelete_mat\OutputVEqvR_YesSean3.mat');

Output=[pt1.Output;pt2.Output;pt3.Output];
Output_dataNames=pt1.Output_dataNames;
save('noDelete_mat\Output_combined_EqvRvary_yesSean.mat','Output',"Output_dataNames")