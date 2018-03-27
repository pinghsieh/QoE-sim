% extract data from figures
clear;
open('Un_Toff=50_80%load.fig');
nRun = 10;
%h = gcf; %current figure handle
%open('example.fig');
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
x = cell2mat(xdata);
y = cell2mat(ydata);
y = (flipud(y))/nRun;
createfigure_all_Un(x(1,:), y);
%axesObjs = get(h, 'Children');  %axes handles
%dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
%ydata = get(dataObjs, 'YData');
