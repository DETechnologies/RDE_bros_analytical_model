%% Init 
close all force
clear 
clc

load("Output_data_stage1.mat")

%% validate that the data is prepared the right way
o=size(Output);
on=size(Output_dataNames);
if o(2) ~= on(2)
    disp("Something wrong with size of data here; plz check")
end 
%% UI dropdown
fig=uifigure("Name","Analytical Model Parameter Selector");
f = uigridlayout(fig,[2 3]);
f.RowHeight = {'3x','1x'};
f.ColumnWidth = {'1x','1x','1x'};

xval=uidropdown(f,"Items",[Output_dataNames(1,:)]);

yval=uidropdown(f,"Items",[Output_dataNames(1,:)]);

zval=uidropdown(f,"Items",[Output_dataNames(1,:)]);

b=uibutton(f,"Text","Plot","ButtonPushedFcn",@(src,event) plot(src,xval,yval,zval,Output_dataNames,Output));
b.Layout.Row = 2;
b.Layout.Column = 2;

%% main plot
%plot these variable : 
function plot(src,xval,yval,zval,Output_dataNames,Output)
    xIndex = find(contains(Output_dataNames(1,:),xval.Value));
    yIndex = find(contains(Output_dataNames(1,:),yval.Value));
    zIndex = find(contains(Output_dataNames(1,:),zval.Value));

    figure("Name","Analytical Model O/P")
    scatter3(Output(:,xIndex),Output(:,yIndex),Output(:,zIndex));
    title("Analytical Model O/P")

    xlabel(xval.Value)
    ylabel(yval.Value)
    zlabel(zval.Value)
end 

