pt1=load('PlottingNice\AnalyticalModelCalculator_March26Results.mat');
pt2=load('PlottingNice\AnalyticalModelCalculator_March26Results_pt2.mat');
pt3=load('PlottingNice\AnalyticalModelCalculator_March26Results_pt3.mat');
pt4=load('PlottingNice\AnalyticalModelCalculator_March26Results_pt4.mat');

Outputs=[pt1.Outputs;pt2.Outputs;pt3.Outputs;pt4.Outputs;];
% Output_dataNames=pt1.Output_dataNames;
save('PlottingNice\combined.mat','Outputs')