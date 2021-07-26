clc
close all
clear all

% syms rxf rxr ry real

rxf = 1;
rxr = 1.6;
ryl = 0.8;
ryr = 1;

A = [1 1 1 1
    -rxf -rxf rxr rxr
    ry -ry -ry ry];

pinvA = pinv(A);

rz = 0.7;

m = 1400;
g = 9.81;

CD = 0*1.2;
rho = 1.225;
S = 1.8*1.4;

R0 = 100;

v = [20:20:100]/3.6; % linspace(10,100)
beta = linspace(-60,20)/180*pi; %[-45:15:0]

c = cell(length(v),length(beta));
N = cell(length(v),length(beta));
for i = 1:length(v)
    for j = 1:length(beta)
        c{i,j} = [m*g;
            (-m*(v(i)^2)/R0*sin(beta(j))+1/2*rho*S*CD*v(i)^2*cos(beta(j)))*rz;
            (-m*(v(i)^2)/R0*cos(beta(j))+1/2*rho*S*CD*v(i)^2*sin(beta(j)))*rz];
        N{i,j} = pinvA*c{i,j};
    end
end


YMatrix1 = [];
YMatrix2 = [];
YMatrix3 = [];
YMatrix4 = [];
for j = 1:length(v)
    temp1 = zeros(1,length(beta)); %length(v)
    temp2 = zeros(1,length(beta));
    temp3 = zeros(1,length(beta));
    temp4 = zeros(1,length(beta));
    for i = 1:length(beta) 
        n = N{j,i};
        temp1(i) = n(1);
        temp2(i) = n(2);
        temp3(i) = n(3);
        temp4(i) = n(4);
    end
    YMatrix1 = [YMatrix1; temp1];
    YMatrix2 = [YMatrix2; temp2];
    YMatrix3 = [YMatrix3; temp3];
    YMatrix4 = [YMatrix4; temp4];
end

max1 = max(max(YMatrix1));
max2 = max(max(YMatrix2));
max3 = max(max(YMatrix3));
max4 = max(max(YMatrix4));
MAX = max([max1, max2, max3, max4]);

% plots
LT = [];
LC = [];
tick = 1;
Nf = length(v);
legenda = cell(Nf,1);
for i = 1:Nf
    %legenda{i} = ['$\beta = $',num2str(beta(i)*180/pi),'$^o$']; 
    legenda{i} = ['$v = $',num2str(v(i)*3.6),'km/h']; 
end

% assi = {'$v$ [km/h]','$N_1$ [N]','Front-Left Tyre Load'};
assi = {'$\beta$ [deg]','$N_1$ [N]','Front-Left Tyre Load'};
createfigure(Nf, beta*180/pi, YMatrix1, assi, legenda, tick, LT, LC)
ylim([0 1.10*MAX])

% assi = {'$v$ [km/h]','$N_2$ [N]','Front-Right Tyre Load'};
assi = {'$\beta$ [deg]','$N_2$ [N]','Front-Right Tyre Load'};
createfigure(Nf, beta*180/pi, YMatrix2, assi, legenda, tick, LT, LC)
ylim([0 1.10*MAX])

%  assi = {'$v$ [km/h]','$N_3$ [N]','Rear-Right Tyre Load'};
assi = {'$\beta$ [deg]','$N_3$ [N]','Rear-Right Tyre Load'};
createfigure(Nf, beta*180/pi, YMatrix3, assi, legenda, tick, LT, LC)
ylim([0 1.10*MAX])

% assi = {'$v$ [km/h]','$N_4$ [N]','Rear-Left Tyre Load'};
assi = {'$\beta$ [deg]','$N_4$ [N]','Rear-Left Tyre Load'};
createfigure(Nf, beta*180/pi, YMatrix4, assi, legenda, tick, LT, LC)
ylim([0 1.10*MAX])
