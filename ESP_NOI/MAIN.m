clc
close all
clear all

%% PARAMETERS

global FzW

g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
Cx0 = -0.3; % [-] drag coefficient 
Cyb = -0.3; % [-] drag coefficient
rk = 1; % dry asphalt
m = 1200; % [kg] vehicle mass
rxf = 1.4;
rxr = -1.6;
ryf = 0.90;
rz =  0.70;
CR = -0.05; % coefficiente attrito resistenza delle ruote
Iz = 1200*(rxr^2+ryf^2)^2;

delta_front = 0*[1 1 0 0]/180*pi; % 3.5 // 10
delta_rear = 0;

rx = [rxf rxf rxr rxr];
ry = [ryf -ryf -ryf ryf];
H = [ones(1,4);
    ry;
    -rx];

%% EQUILIBRIUM POINT

x0 = 0;
y0 = 0;
v0 = 130/3.6; % [m/s] vehicle speed
beta0 = 0;
w0 = 0;
psi0 = 0/180*pi;
theta = 15/180*pi;

pinvH = pinv(H);
FzW = pinvH*(rz/2*rho*S*v0^2*[0; Cyb*sin(beta0); -Cx0*cos(beta0)]+[m*g*cos(theta); -rz*m*v0*w0*cos(beta0); -rz*m*g*sin(theta)-rz*m*v0*w0*sin(beta0)]); % forza ground su routa

lambda0s = fsolve(@long_eq,0);
lambda0 = [lambda0s; lambda0s; 0; 0];

Cf = (FzW(1)+FzW(2))*Partial_mu_long(1,0); %cornering stiffness front     partial_mu_long derivata parziale rispetto a beta
Cr = (FzW(3)+FzW(4))*Partial_mu_long(1,0); % rear

A = [-(Cf+Cr)/(m*v0) -v0-(Cf*rxf+Cr*rxr)/(m*v0)
    -(Cf*rxf+Cr*rxr)/(Iz*v0) -(Cf*rxf^2+Cr*rxr^2)/(Iz*v0)];
B1_delta = Cf*[1/m;
      rxf/Iz];
B1_lambda = [1/m 1/Iz]*[zeros(1,4) Cr;
     [FzW(1)*(-Partial_mu_long(1,lambda0s))*ryf FzW(2)*Partial_mu_long(1,lambda0s)*ryf FzW(3)*Partial_mu_long(1,0)*ryf -FzW(4)*Partial_mu_long(1,0)*ryf Cr*rxr]];

B2eq = [0; 1/Iz];

C1 = [0 1]; % [-] measured output
y0 = w0; % [rad/s] linearisation output

s = (Cf*rxf+Cr*rxr)/(Cf+Cr);
if s > 0
    disp('OVERSTEERED VEHICLE')
    disp(['s = ',num2str(s),' m'])
elseif s < 0
    disp('UNDERSTEERED VEHICLE')
    disp(['s = ',num2str(s),' m'])
else
    disp('NEUTRAL STEERED VEHICLE')
end

%% MU LONGITUDINAL

function mu = mu_long(i,lambda) %12.1d

% 1) Dry asphalt
% 2) wet asphalt
% 3) Snow
% 4) Ice
% 5) Dry Cobblestone
% 6) wet cobblestone
theta1 = [1.28 0.86 0.19 0.05 1.37 0.4];
theta2 = [23.99 33.82 94.13 306.39 6.46 33.71];
theta3 = [0.52 0.35 0.05 0 0.67 0.12];

mu = sign(lambda).*theta1(i).*(1-exp(-abs(lambda)*theta2(i)))-theta3(i)*lambda;
end

%% LONGITUDINAL EQUATION

function f = long_eq(x) % LONGITUDINAL EQUATION

global FzW

rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
Cx0 = -0.3; % [-] drag coefficient 
v0 = 130/3.6; % initial speed
CR = -0.05; 
kind = 1; % kind of road surface

f = 1/2*rho*S*v0^2*Cx0 + FzW.'*[mu_long(kind,x)+CR
                        mu_long(kind,x)+CR
                        CR
                        CR];
                    
% eq. che porta il punto di equilibrio con forza sulle ruote = drag
end

%% PARTIAL DERIVATIVE MU LONGITUDINAL

function mu = Partial_mu_long(i,lambda) % partial_mu_long derivata parziale rispetto a beta (slittamento laterale)

% 1) Dry asphalt
% 2) wet asphalt
% 3) Snow
% 4) Ice
% 5) Dry Cobblestone
% 6) wet cobblestone
theta1 = [1.28 0.86 0.19 0.05 1.37 0.4];
theta2 = [23.99 33.82 94.13 306.39 6.46 33.71];
theta3 = [0.52 0.35 0.05 0 0.67 0.12];

mu = theta1(i).*theta2(i).*exp(-abs(lambda).*theta2(i))-theta3(i);
end

