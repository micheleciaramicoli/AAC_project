function dx = ESPNLsim(t,x)

%g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
CD = 2.2; % [-] drag coefficient

%rk = 1; % dry asphalt
a = 1.5; % [m]
b = 4-a; % [m]
M = 1200; % [kg] vehicle mass
%N = 2*M*g/4; % [N] vertical force
%r = 0.28; % [m] wheel radius
Iz =M/4*(2*a^2+2*b^2); % [kg*m^2] z-axis vehicle inertia

V = x(1);
beta = x(2);
wz = x(3);

Fw1b = [cos(delta) -sin(delta)
        sin(delta)  cos(delta)]*Fw1;
Fw2b = Fw2;
    
Fb = Fw1b + Fw2b;
Mb = Fw1b(2)*a - Fw2b(2)*b;

e = [beta - beta0;
     wz - wz0];
tau = -K*e;

dV = 1/M*[cos(beta) sin(beta)]*Fb-1/M*1/2*rho*S*V^2*CD;
dbeta = 1/(M*V)*[-sin(beta) cos(beta)]*Fb - wz;
dwz = 1/Iz*Mb+1/Iz*tau;

dx = [dV
      dbeta
      dwz];