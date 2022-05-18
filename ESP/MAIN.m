% MAIN.mat

clc 
clear 
close all

%% INPUT VARIABLES

global rho S v0 Cx N_wheels CR
global g Cy rk r_x r_y r_z m J

g   = 9.81;                         % [m/s^2] gravity acceleration
rho = 1.225;                        % [kg/m^3] air density
Cx  = 0.3;                          % drag coefficient 
Cy  = 0.3;                          % drag coefficient
rk  = 3;                            % dry asphalt
v0  = 120/3.6;

r_x = [1.6 1.6 -1.4 -1.4];
r_y = [0.9 -0.9 0.9 -0.9];
r_z = 0.70;
S   = 2.4;                          % [m^2] cross surface
m   = 1200;                         % [kg] vehicle mass
J   = m*(r_x(4)^2+r_y(1)^2)^2;      % Vehicle inertia

CR = 0.05;                          % coefficiente attrito resistenza delle ruote

%% parameterns at eq point

omega           = 0;
beta_body       = 0;
v_x_body        = v0 * cos(beta_body);    % speed of the body in the x direction
v_y_body        = v0 * sin(beta_body);    % speed of the body in the y direction
delta_w0        = [0; 0; 0; 0];                 % steering angles directly of the wheels (non c'è il volante)
delta_f         = delta_w0(1);           
delta_r         = delta_w0(3);

N               = m*g;

H               = [ones(1,4);
                    -r_x;
                     r_y   ];

N_wheels        = [3320.6;3320.6;2364.8;2364.8];

% N_wheels        = pinv(H) * [          N
%                               0.5*rho*(v0^2)*S*Cx*r_z
%                                        0               ];
%                                   

u_delta_f0      = 50*pi/180;

lambda0         = fsolve(@long_eq,0);
lambdas0_w      = [lambda0; lambda0; 0; 0];

uFB0            = [delta_r; lambdas0_w];

x0              = [0; 0; v0; 0; 0; 0];

Cf          = (N_wheels(1) + N_wheels(2)) * Partial_mu_long(rk,0);
Cr          = (N_wheels(3) + N_wheels(4)) * Partial_mu_long(rk,0);

eta         =  - N * ((Cf * r_x(1) + Cr * r_x(3))/(Cf * Cr *(r_x(1) - r_x(3))))  % r_x(3) è negativo

if eta < 0
V_critical  = (sqrt((-g/eta)*(r_x(1) - r_x(3)))) * 3.6  
end

%% FUNCTIONS

function f = long_eq(x) % LONGITUDINAL EQUATION

global rho S v0 Cx N_wheels CR

kind = 1;

f = -1/2*rho*S*v0^2*Cx + N_wheels.'*[mu_long(kind,x)-CR
                        mu_long(kind,x)-CR
                        -CR
                        -CR];
                    
% eq. che porta il punto di equilibrio con forza sulle ruote = drag
end

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