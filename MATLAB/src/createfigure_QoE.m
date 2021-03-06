function createfigure_QoE(X1, D1, D2, D3, D4, D5)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 09-Jun-2015 11:00:18

% Create figure
figure1 = figure('Color',[1 1 1]);
YMatrix1 = vertcat(D1, D2, D3, D4, D5);
% Create axes
axes1 = axes('Parent',figure1,'LineWidth',4,'FontWeight','bold',...
    'FontSize',20);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'Marker','o');
set(plot1(1),'Marker','*','Color',[1 0 0]);
set(plot1(2),'Color',[0 1 0]);
set(plot1(3),'Color',[0 0 1]);
set(plot1(4),'Color',[1 0 1]);
set(plot1(5),'Marker','*','Color',[0 1 1]);

% Create xlabel
xlabel({'Time (# of slots)'},'FontWeight','bold','FontSize',20);

% Create ylabel
ylabel({'D_1(t)'},'FontWeight','bold','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.156109234234235 0.58921036769138 0.105574324324324 0.30379746835443],...
    'LineWidth',4);

