clc
close all
clear all

% A matrix ESP
V = 90/3.6;
b = 2/180*pi;
R = 75;
rk = 1;
delta1 = 10/180*pi;
delta2 = 10/180*pi;
wz = V/R;


rho = 1.225; % [kg/m^3] air density
g = 9.81; % [m/s^2] gravity acceleration

J = 1; % [kg m^2] wheel inertia
m = 1200; % [kg] vehicle mass
r = 0.28; % [m] wheel radius
CX0 = -0.34; % [-] drag coeffcient
CYb = -0.034; % [-] drag coeffcient
S = 2.4; % [m^2] vehicle cross section
CR = 0.015; % [-] rolling resistance coefficient
Iz = 4*(2^2)*m/4;


n = 4;
ry = 1.2;
rx1 = 1.5;
rx3 = 4-rx1;
rz = 0.70;
ryV = [ry;-ry;-ry;ry];
rxV = [rx1; rx1; -rx3; -rx3];

w1 = wz*(-ryV(1))+V*cos(b);
w2 = wz*(-ryV(2))+V*cos(b);
w3 = wz*(-ryV(3))+V*cos(b);
w4 = wz*(-ryV(4))+V*cos(b);

deltaV = [delta1; delta2; 0 ;0];
wV = [w1, w2, w3, w4];


MM = diag([m*V, Iz]);

H = [ones(1,4);
    ryV';
    -rxV' ];

FzB =  (H')*inv(H*H')*(rz/2*rho*S*V^2*[0; CX0; -CYb*b]+[m*g; -rz*m*V^2/R*cos(b); -rz*m*V^2/R*sin(b)]);

FaB = 1/2*rho*V^2*S*[CX0; CYb*b; 0];

FB = FaB;
for i = 1:n
    
    FwzW = FzB(i);
    ry = ryV(i);
    rx = rxV(i);
    d = deltaV(i);
    w = wV(i);
    
    Rdt = [cos(d) sin(d);
        -sin(d) cos(d)];
    VW = Rdt*(wz*[-ry;rx]+V*[sin(b); cos(b)]);
    VxW = VW (1);
    VyW = VW (2);
    vmax = sqrt(VyW^2+(max(w*r,VxW))^2);
    LL = (w*r-VxW)/vmax;
    LS = -VyW/vmax;
    LT = sqrt(LL^2+LS^2);
    mu = mu_long(rk,LT);
    muL = LL/LT*mu;
    muS = LS/LT*mu;
    FwW = FwzW*[muL; muS; 1];
    Rd = [cos(d) -sin(d) 0
        sin(d) cos(d) 0
        0 0 1];
    FwB = Rd*FwW;
    FB = FB+FwB;
end

FxB = FB(1);
FyB = FB(2);
FzB = FB(3);

% A11 = cos(b)*(dFyBdb-FxB)-sin(b)*(FyB+dFxBdb);
% A12 = -m*V0-sin(b)*dFyBdwz-sin(b)*dFxBdwz;
% A21 = dtdb;
% A22 = dtdwz;
% 
% A = MM\[A11 A12;
%     A21 A22];