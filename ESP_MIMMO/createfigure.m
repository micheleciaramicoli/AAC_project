function createfigure(N, XX, YMatrix, assi, legenda, tick, LT, LC)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 31-Jan-2020 15:32:06

LW = 2;
TITLE_SIZE = 30;
TEXT_SIZE = 24;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

[a,~] = size(XX);

% Create multiple lines using matrix input to plot
for i = 1:N
    if a == 1
        X = XX;
    else
        X = XX(i,:);
    end
    if isempty(LT) == 1
        if isempty(LC) ==1
            plot1(i) = plot(X,YMatrix(i,:),'LineWidth',LW);
        else
            plot1(i) = plot(X,YMatrix(i,:),'LineWidth',LW,'color',LC{i});
        end
    else
        if isempty(LC) ==1
            plot1(i) = plot(X,YMatrix(i,:),LT{i},'LineWidth',LW);
        else
            plot1(i) = plot(X,YMatrix(i,:),LT{i},'LineWidth',LW,'color',LC{i});
        end
    end
    if isempty(legenda) == 1
    else
        set(plot1(i),'DisplayName',legenda{i});
    end
end

% Create ylabel
ylabel(assi{2},'Interpreter','latex', 'FontSize',TEXT_SIZE);

% Create xlabel
xlabel(assi{1},'Interpreter','latex', 'FontSize',TEXT_SIZE);


box(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',TEXT_SIZE);
% Create legend
if isempty(legenda) == 1
else
    legend1 = legend(axes1,'show');
    set(legend1,'Interpreter','latex');
end

if tick == 0
    set(gca,'YTick',0)
    set(gca,'XTick',0)
    
end
grid(axes1,'on');

% Create title
title(assi{3},'FontWeight','bold','FontSize',TITLE_SIZE,'Interpreter','latex');