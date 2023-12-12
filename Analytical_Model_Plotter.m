%% Init 
close all force
clear 
clc

load("Output_data_dec10_new.mat")

%% validate that the data is prepared the right way
o=size(Output);
on=size(Output_dataNames);
if o(2) ~= on(2)
    disp("Something wrong with size of data here; plz check")
end 
%% UI dropdown - 3d plot
fig1=uifigure("Name","Analytical Model 3D Plotter");
f1 = uigridlayout(fig1,[2 3]);
f1.RowHeight = {'3x','2x'};
f1.ColumnWidth = {'1x','1x','1x'};

xval=uidropdown(f1,"Items",[Output_dataNames(1,:)]);

yval=uidropdown(f1,"Items",[Output_dataNames(1,:)]);

zval=uidropdown(f1,"Items",[Output_dataNames(1,:)]);    

b=uibutton(f1,"Text","Plot","ButtonPushedFcn",@(src,event) plot3D(src,xval,yval,zval,Output_dataNames,Output));
b.Layout.Row = 2;
b.Layout.Column = 2;

%% UI dropdown - 2D plot
fig2=uifigure("Name","Analytical Model 2D Plotter");
f2 = uigridlayout(fig2,[3 2]);
f3.RowHeight = {'3x','2x'};
f3.ColumnWidth = {'1x','1x'};

xval=uidropdown(f2,"Items",[Output_dataNames(1,:)]);

yval=uidropdown(f2,"Items",[Output_dataNames(1,:)]);

cb=uicheckbox(f2,"Text","Two y-axis?");
yval2=uidropdown(f2,"Items",[Output_dataNames(1,:)]);

b=uibutton(f2,"Text","Plot","ButtonPushedFcn",@(src,event) plot2D(src,xval,yval,yval2,Output_dataNames,Output,cb));
b.Layout.Row = 3;
b.Layout.Column = 2;



%% plot functions;
%plot these variable : 
function plot2D(src,xval,yval,yyvalueRight,Output_dataNames,Output,secondAx)
    xIndex = find(contains(Output_dataNames(1,:),xval.Value));
    yIndex = find(contains(Output_dataNames(1,:),yval.Value));

    figure("Name","Analytical Model 2D Plot O/P")
    scatter(Output(:,xIndex),Output(:,yIndex));
    title("Analytical Model O/P")

    xlabel(xval.Value)
    ylabel(yval.Value)

    if secondAx.Value == 1
        yyindexRight = find(contains(Output_dataNames(1,:),yyvalueRight.Value));
        yyright=Output(:,yyindexRight);
        xxval=Output(:,xIndex);
        addAxis2Plot(xxval,yyright);
        ylabel(yyvalueRight.Value)
    end 
 
end 

function addAxis2Plot(xxval,yyright)
    yyaxis right 
    scatter(xxval,yyright);
end 


%plot these variable : 
function plot3D(src,xval,yval,zval,Output_dataNames,Output)
    xIndex = find(contains(Output_dataNames(1,:),xval.Value));
    yIndex = find(contains(Output_dataNames(1,:),yval.Value));
    zIndex = find(contains(Output_dataNames(1,:),zval.Value));

    figure("Name","Analytical Model 3D Plot O/P")
    scatter3(Output(:,xIndex),Output(:,yIndex),Output(:,zIndex));
    title("Analytical Model O/P")

    xlabel(xval.Value)
    ylabel(yval.Value)
    zlabel(zval.Value)
end 


