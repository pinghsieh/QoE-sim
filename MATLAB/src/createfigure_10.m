function createfigure_10(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 05-Nov-2017 21:08:40

% Create figure
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','1','Color',[1 0 0]);
set(plot1(2),'DisplayName','2','Color',[0 1 0]);
set(plot1(3),'DisplayName','3','Color',[0 0 1]);
set(plot1(4),'DisplayName','4','Color',[1 0 1]);
set(plot1(5),'DisplayName','5','Color',[1 1 0]);
set(plot1(6),'DisplayName','6','Color',[0 1 1]);
set(plot1(7),'DisplayName','7','Color',[1 0 0.3]);
set(plot1(8),'DisplayName','8','Color',[0.5 0.8 0]);
set(plot1(9),'DisplayName','9','Color',[0.3 0 1]);
set(plot1(10),'DisplayName','10','Color',[0.8 0.2 0]);

% Create xlabel
xlabel({'Time (slot)'});

% Create ylabel
ylabel({'D_n(t)'});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[4326.03686635945 6468.89400921659]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[6.69072164948454 8.25773195876289]);
box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',17,'LineWidth',2);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.152232142857143 0.447619047619048 0.104464285714286 0.460714285714286]);

