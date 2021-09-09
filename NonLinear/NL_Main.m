clc
close all
clear all

%% PARAMETERS

g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
Cx0 = -0.3; % [-] drag coefficient (IN THE NOTES ONLY Cd ON X AXIS)
Cyb = -0.3; % [-] drag coefficient
rk = 1; % dry asphalt
m = 1200; % [kg] vehicle mass
rxf = 1.4;
rxr = -0.4;
ryf = 0.90;
rz =  0.70;
CR = -0.05;
Iz = 1200*(rxr^2 + ryf^2)^2;

rx = [rxf rxf rxr rxr];
ry = [ryf -ryf -ryf ryf];
H = [ones(1,4); %3x4
    ry;
    -rx];

%% EQUILIBRIUM POINT (linearization trajectory)

x0 = 0;
y0 = 0;
v0 = 130/3.6;           % [m/s] vehicle speed from GPS independently from the direction
beta0 = 0;              % side-slip angle, 0 for the linearization trajectory (angolo tra vettore velocità e asse x riferito al body)
v0x = v0*cos(beta0);    % actual speed in x from the measured GPS speed
w0 = 0;                 % yaw angular speed
psi0 = 0/180*pi;        % yaw angle, 0 for the linearization trajactory
theta = 0;              % 15/180*pi road steepness
delta_r0 = 0;           % equilibrium rear wheels = non steered

% questo per trovare le forze verticali sulle 4 ruote (N1234);
% calcolate sempre rispetto ai valori 
pinvH = pinv(H); %pseudoinverse will be 4x3
FzW = pinvH*(rz/2*rho*S*v0^2*[0; Cyb*sin(beta0); -Cx0*cos(beta0)] + [m*g*cos(theta); -rz*m*v0*w0*cos(beta0); -rz*m*g*sin(theta)-rz*m*v0*w0*sin(beta0)]);

REF = [rho, S, v0, CR, Cx0];
f = @(x)long_eq(x,REF,FzW); 

lambda0s = fsolve(f,0); %0.0061 @130, 0.0208 @270
lambda0 = [lambda0s; lambda0s; 0; 0]; 

u0 = [lambda0; delta_r0]; 

Cf = (FzW(1)+FzW(2))*Partial_mu_long(1,0);
Cr = (FzW(3)+FzW(4))*Partial_mu_long(1,0);
s = (Cf*rxf + Cr*rxr)/(Cf + Cr);

eta = -m*g*(Cf*rxf + Cr*rxr)/(Cf*Cr*(rxf-rxr));

crit_Speed = sqrt(g*(rxf-rxr)/abs(eta));

if s > 0
    disp('OVERSTEERED VEHICLE')
    disp(['s = ',num2str(s),' m'])
elseif s < 0
    disp('UNDERSTEERED VEHICLE')
    disp(['s = ',num2str(s),' m'])
else
    disp('NEUTRAL STEERED VEHICLE')
end

disp(['eta = ',num2str(eta)])
disp(['Critical Speed = ',num2str(crit_Speed*3.6),' km/h'])


%% REACHABILITY

% B_eq = [0 Cr/m
%     1/Iz Cr*(rxf/Iz)];
% R = [B_eq A*B_eq];
% 
% if rank(R) == 2
%     disp('FULLY REACHABLE')
% else
%     disp('NOT FULLY REACHABLE')
% end
% 
% %% OBSERVABILITY
% 
% O = [C1
%     C1*A];
% 
% if rank(O) == 2
%     disp('FULLY OBSERVABLE')
% else
%     disp('NOT FULLY OBSERVABLE')
% end

% %% CONTROL
% 
% Kr = [1 1
%     1 1];
% 
% de_tau = [B1_Ctrl(2,1) B1_Ctrl(2,2) B1_Ctrl(2,3) B1_Ctrl(2,4)];
% pinv_de_tau = pinv(de_tau);
% 

%% LONGITUDINAL EQUILIBRIUM EQUATION

function f = long_eq(x,REF,FzW)

kind = 1;

f = 1/2*REF(1)*REF(2)*REF(3)^2*REF(5)+(FzW.')*[mu_long(kind,x)+REF(4)
                                             mu_long(kind,x)+REF(4)
                                             REF(4)
                                             REF(4)];
end

%% MU

function mu = mu_long(i,lambda)

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

%% Partial derivative of mu

function mu = Partial_mu_long(i,lambda)

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

%% Steering Wheel Function

function delta_f = SteeringWheel(t)

ti = 5; % [s] inital time
tf = 7; % [s] final time
A = 5/180*pi; %[rad] steering wheel maximum amplitude
%FaB = 1/2*rho*S*v^2*Cd;

type = 0;
% 0 = sin
% 1 = step

switch type
    case 0
        % SIN COMPUTATION
        DT = tf-ti; % [s] steering action duration
        w = 2*pi/DT;
        
        if t >= ti && t <= tf
            delta_f = A*sin(w*(t-ti));
        else
            delta_f = 0;
        end
    case 1
        if t >= ti 
            delta_f = A;
        else
            delta_f = 0;
        end
end
end
