% MAIN.mat

clc 
clear 
close all

%% VARIABLES INPUT

global g rho Cx Cy rk v
global r_x r_y r_z S m J CR

g   = 9.81;                   % [m/s^2] gravity acceleration
rho = 1.225;                  % [kg/m^3] air density
Cx  = 0.3;                    % drag coefficient 
Cy  = 0.3;                    % drag coefficient
rk  = 1;                      % dry asphalt
v   = 50/3.6;

r_x = [1.4 1.4 -1.6 -1.6];
r_y = [0.9 -0.9 0.9 -0.9];
r_z =  0.70;
S   = 2.4;                          % [m^2] cross surface
m   = 1200;                         % [kg] vehicle mass
J   = m*(r_x(4)^2+r_y(1)^2)^2;      % Vehicle inertia

CR = -0.05;                 % coefficiente attrito resistenza delle ruote

%% EQUATIONS

omega           = 0;
beta_body       = 0;
v_x_body        = v * cos(beta_body);    % speed of the body in the x direction
v_y_body        = v * sin(beta_body);    % speed of the body in the y direction
delta_w0        = [0 0 0 0];        % steering angles

N               = m*g;

H               = [ ones(1,4)
                    -r_x
                     r_y   ];

N_wheels        = pinv(H) * [        N
                               0.5*rho*(v^2)*S*Cx*r_z
                                      0                ];

lambdas0        = ;


v_x_wheels      = [1 1 1 1];        % speed of the wheels        
v_y_wheels      = [1 1 1 1];        % speed of the wheels    
tau             = [1 1 1 1];        % torques
f_x             = [1 1 1 1];
f_y             = [1 1 1 1];
f_x_wheels      = [1 1 1 1];
f_y_wheels      = [1 1 1 1];
lambdas         = [1 1 1 1];        % longitudinal slip ratios
betas_wheels    = [1 1 1 1];        % side slip angles

for i = 1:4
    Rot_wtob        = [+cos(delta_w(i)) sin(delta_w(i))
                       -sin(delta_w(i)) cos(delta_w(i))];
    
    V_wheels        = Rot_wtob*([-omega*r_y(i) ; omega*r_x(i)] + [v_x_body ; v_y_body]);
    v_x_wheels(i)   = V_wheels(1);
    v_y_wheels(i)   = V_wheels(2);

    betas_wheels(i) = atan(v_y_wheels(i)/v_x_wheels(i));

    V_wheels        = [0;0];         % Reset value of V_wheels
end

for i = 1:4
    mu_lambda   = mu_long(rk, lambdas(i));              % Dry asphalt
    mu_beta     = mu_long(rk, betas_wheels(i));         % Dry asphalt
end    

for i = 1:4
    f_x_wheels(i) = + N_wheels(i)*mu_lambda;
    f_y_wheels(i) = - N_wheels(i)*mu_beta;
end 

for i = 1:4
    f_x     = f_x + (cos(delta_w(i))*f_x_wheels(i) -sin(delta_w(i))*f_y_wheels(i));
    f_y     = f_y + (sin(delta_w(i))*f_x_wheels(i) +cos(delta_w(i))*f_y_wheels(i));
    tau     = tau + (sin(delta_w(i))*f_y_wheels(i) -cos(delta_w(i))*f_x_wheels(i))*r_y(i) + (sin(delta_w(i))*f_x_wheels(i) +cos(delta_w(i))*f_y_wheels(i))*r_x(i);
end 

v_xdot      = (f_x + (m*omega*v_y_body) - (0.5*rho*S*(v^2)*Cx))/m;
v_ydot      = (f_y - (m*omega*v_x_body))/m;
omega_dot   = tau/J;

x_dot       = [v_xdot
               v_ydot 
               omega_dot];      % Derivata degli stati


















    
