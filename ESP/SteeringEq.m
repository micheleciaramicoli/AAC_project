function f = SteeringEq(x)

global rk
global V delta a b N

g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
CD = 2.2; % [-] drag coefficient

%V = 90/3.6; % [m/s] vehicle speed
%a = 1.5; % [m]
%b = 4-a; % [m]
%delta = 0/180*pi; % [rad] mean steerign angle
M = 1200; % [kg] vehicle mass
% N = M*g/4; % [N] vertical force
r = 0.28; % [m] wheel radius

R = x(1);
%tau = x(2);
beta = x(2);
w1 = x(3);
%betaw1 = % x(5);

% R1 = sin(pi/2+beta)/sin(pi/2-betaw1-delta)*R;
R1 = sqrt(a^2 + R^2 -2*a*R*cos(pi/2+beta));
betaw1 = (-delta-asin(R/R1*sin(pi/2+beta))+pi/2); % x(6);
% R2 = sin(pi/2-beta)/sin(pi/2-betaw2)*R;
R2 = sqrt(b^2 + R^2 -2*b*R*cos(pi/2-beta));
betaw2 = (asin(R/R2*sin(pi/2-beta))-pi/2); % x(6);

Vw1 = V/R*R1;
Vw2 = V/R*R2;

lambda1L = (Vw1*cos(betaw1)-w1*r)/max(abs(w1*r),Vw1*cos(betaw1));
lambda1S = (Vw1*sin(betaw1))/max(abs(w1*r),Vw1*cos(betaw1));
lambda1TOT= sqrt(lambda1L^2+lambda1S^2);

w2 = 1/r*Vw2*cos(betaw2);
lambda2L = 0; %(Vw2*cos(betaw2)-w2*r)/(w2*r);
lambda2S = (Vw2*sin(betaw2))/(w2*r);
lambda2TOT= sqrt(lambda2L^2+lambda2S^2);

if lambda1TOT == 0;
    Fxw1 = 0;
    Fyw1 = 0;
else
    Fxw1 = -lambda1L/lambda1TOT*r*N*mu_long(rk,lambda1TOT);
    Fyw1 = -lambda1S/lambda1TOT*r*N*mu_long(rk,lambda1TOT);
end
if lambda2TOT == 0;
    Fyw2 = 0;
else
    Fyw2 = -lambda2S/lambda2TOT*r*N*mu_long(rk,lambda2TOT);
end

Fxb = Fxw1*cos(delta)-Fyw1*sin(delta);
Fyb = Fxw1*sin(delta)+Fyw1*cos(delta)+Fyw2;

f = [ %sin(pi/2+betaw2)/R - sin(-betaw2+beta)/b;% lower triangle
    %sin(pi/2-betaw1-delta)/R - sin(betaw1+delta+beta)/a;% upper triangle
    M*V^2/R + (sin(beta)*Fxb-cos(beta)*Fyb); % ayv
    1/2*rho*S*V^2*CD-cos(beta)*Fxb-sin(beta)*Fyb; % axv
    a*(Fxw1*sin(delta)+Fyw1*cos(delta))-Fyw2*b; % z-momentum
    %-Fxw1/2+r*tau
    ];