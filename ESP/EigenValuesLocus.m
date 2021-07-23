clc
close all
clear all

Cf = 26356*2;
m = 917;
g = 9.81;
Iz = 1128;
rxf = 0.91;
rxr = -1.64;
Cr = 28648*2; %  -Cf*rxf/rxr; %

inversion = 1;

if inversion == 1
    temp=Cf;
    Cf=Cr;
    Cr=temp;
    temp=rxf;
    rxf=-rxr;
    rxr=-temp;
end

v0 = linspace(10,150,1000)/3.6;

eta = -m*g*(Cf*rxf+Cr*rxr)/(Cr*Cf*(rxf-rxr));
if eta > 1e-10
    disp('OVERSTEERED VEHICLE')
    titleS = '>0$';
elseif eta < -1e-10
    disp('UNDERSTEERED VEHICLE')
    vcr = sqrt(-g/eta*(rxf-rxr));
    titleS = '<0$';
else
    disp('NEUTRAL STEERED VEHICLE')
    titleS = '=0$';
end
N = length(v0);
for i = 1:N
A11 = -(Cf+Cr)/(m*v0(i));
A12 = -v0(i)-(Cf*rxf+Cr*rxr)/(m*v0(i));
A21 = -(Cf*rxf+Cr*rxr)/(Iz*v0(i));
A22 = -(Cf*rxf^2+Cr*rxr^2)/(Iz*v0(i));
temp = eig([A11 A12
            A21 A22]);        
e1r(1,i) = real(temp(1));
e1i(1,i) = imag(temp(1));
e2r(1,i) = real(temp(2));
e2i(1,i) = imag(temp(2));
end

% plots
LT = [];
LC = [];
tick = 0;
Nf = 2;
assi = {'$v_0$ ','Real(s)',['Eigenvalues Locus for $\eta',titleS]};
legenda = {'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
YMatrix = [e1r; e2r];
createfigure(Nf, v0*3.6.', YMatrix, assi, legenda, tick, LT, LC)
grid on
Nf = 2;
assi = {'$v_0$ ','Imag(s)',['Eigenvalues Locus for $\eta',titleS]};
legenda = {'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
YMatrix = [e1i; e2i];
createfigure(Nf, v0*3.6.', YMatrix, assi, legenda, tick, LT, LC)


%% LINEAR LYAPUNOV TEST
v0 = 100/3.6;
A11 = -(Cf+Cr)/(m*v0);
A12 = -v0-(Cf*rxf+Cr*rxr)/(m*v0);
A21 = -(Cf*rxf+Cr*rxr)/(Iz*v0);
A22 = -(Cf*rxf^2+Cr*rxr^2)/(Iz*v0);
A = [A11 A12; A21 A22];
A = [A zeros(2,1)
     0 1 0];
 
B = [0;0;0];
C = eye(3);
D = zeros(3,1);
sys = ss(A,B,C,D);
t = linspace(0,1,1000);
u = zeros(1,length(t));
x0 = [v0; 0; 0];
[y,t,x] = lsim(sys,u,t,x0);
p(:,1) = zeros(2,1);
for i = 1:length(t)-1
    p(:,i+1) = p(:,i) + (t(i+1)-t(i))*[cos(x(i,3)) -sin(x(i,3)); sin(x(i,3))  cos(x(i,3))]*[v0; x(i,1)];
end

Nf = 2;
assi = {'$t$ ','$\tilde{\bf x}$',['Time History for $\eta',titleS]};
legenda = {'$\tilde{v}_y$','$\tilde{\omega}$'}; %{'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
YMatrix = x(:,1:2).';
createfigure(Nf, t.', YMatrix, assi, legenda, tick, LT, LC)

% Nf = 1;
% assi = {'$t$ ','$\tilde{\omega}$',['Time History for $\eta',titleS]};
% legenda = []; %{'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
% YMatrix = x(:,2).';
% createfigure(Nf, t, YMatrix, assi, legenda, tick, LT, LC)

Nf = 1;
assi = {'$t$ ','$\tilde{\psi}$',['Time History for $\eta',titleS]};
legenda = []; %{'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
YMatrix = x(:,3).';
createfigure(Nf, t.', YMatrix, assi, legenda, tick, LT, LC)

Nf = 1;
assi = {'${p}_x$','${p}_y$',['Time History for $\eta',titleS]};
legenda = []; %{'$s_1$','$s_2$'}; %{['CBP = ',num2str(CBPv(1))],['CBP = ',num2str(CBPv(2))],['CBP = ',num2str(CBPv(3))],['CBP = ',num2str(CBPv(4))],['CBP = ',num2str(CBPv(5))],['CBP = ',num2str(CBPv(6))],'Closed Loop'};
YMatrix = p(2,:);
createfigure(Nf, p(1,:), YMatrix, assi, legenda, tick, LT, LC)