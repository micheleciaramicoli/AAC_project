clc
close all
clear all

%% PARAMETERS

global m rho S g rxf rxr ryf rz Iz theta
global Cx0 Cyb CR rk pinvH
global lambda0
global tmin tmax
global Klqr CL ICE A B1 Ce y0

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
CR = -0.05;
Iz = 1200*(rxr^2+ryf^2)^2;

delta = 0*[1 1 0 0]/180*pi; % 3.5 // 10

rx = [rxf rxf rxr rxr];
ry = [ryf -ryf -ryf ryf];
H = [ones(1,4);
    ry;
    -rx];

%% EQUILIBRIUM POINT

global v0 FzW

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
A11 = -(Cf+Cr)/(m*v0);
A12 = -v0-(Cf*rxf+Cr*rxr)/(m*v0);
A21 = -(Cf*rxf+Cr*rxr)/(Iz*v0);
A22 = -(Cf*rxf^2+Cr*rxr^2)/(Iz*v0);
A = [A11 A12
    A21 A22];
B1 = Cf*[1/m;
      rxf/Iz];
B2 = ryf/Iz*[zeros(1,4);
             FzW.*[-Partial_mu_long(1,lambda0s) Partial_mu_long(1,lambda0s) Partial_mu_long(1,0) -Partial_mu_long(1,0)]]; %this is B1lambda
         
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

%% STABILISER
% CONTROL
% Resonanace condition check
Ce = [0 1];
Mcheck = [A B2eq
          Ce 0];
if rank(Mcheck) == 3
    disp('RESONANCE CONDITION VERIFIED')
else
    disp('ERROR: RESONANCE CONDITION NOT VERIFIED')
end
Ac = [Ce 0;
      A zeros(2,1)];
Bc = [0;
      B2eq];
 
Qx = inv(diag([(2)^2 (20/180*pi)^2]));
Qeta = 0.1*Ce*Qx*Ce.';
Qp = 1/3*blkdiag(Qeta,Qx);
Req = inv((ryf*2500)^2);
[Klqr,Sinfty,eigACL] = lqr(Ac,Bc,Qp,Req);
disp('- LQR GAIN = ');disp(-Klqr);

% OBSERVER
sys = ss(A, [B1 B2eq eye(2)], C1, [zeros(1,2) zeros(1,2)]);
% process noise
Qd11 = diag([0.1, 0.1]);
% sensor noise
std_gyro = 1/3*1/180*pi; % [m/s] standard deviation of the gyroscope
Rd = std_gyro^2; % [(m/s)^2] variance of gyroscope 
% cross terms
Qd12 = zeros(2,1);
global kest
kest = kalman(sys,Qd11,Rd,Qd12);

%% OPEN LOOP SIMULATION
        
tspan = [0:0.1:30];
tmin = max(tspan)/3 -1;
tmax = max(tspan)/3 +1;
%options = odeset('MaxStep',0.01);
ICE = 0;
CL = 0;
[t,y] = ode23(@ESPNLSim_v1,tspan,[x0; y0; v0; beta0+0*45/180*pi; w0; psi0; zeros(2,1); zeros(2,1)]);
CL = 1;
[tCL,yCL] = ode23(@ESPNLSim_v1,tspan,[x0; y0; v0; beta0+0*45/180*pi; w0; psi0; zeros(2,1);zeros(2,1)]);

% ICE = 1;
% [tICE,yICE] = ode23(@ESPNLSim_v1,tspan,[x0; y0; v0; beta0; w0; psi0],options);

w_ref = zeros(1,length(tCL));
delta_f = w_ref;
tau_eq = w_ref;
for i = 1:length(tCL)
    delta = SteeringWheel(tCL(i));
    delta_f(i) = delta(1); 
    xref = ReferenceGenerator(delta);
    w_ref(i) = xref(2);
    tau_eq(i) = -Klqr(2:3)*([yCL(i,3)*sin(yCL(i,4));yCL(i,5)]-xref)-Klqr(1)*Ce*yCL(i,7:8).';
end

%% plots
LT = [];
LC = [];
tick = 1;
Nf = 2;
assi = {'$t$ [s]','$\tilde{v}_y$ [m/s]','Side speed'};
legenda = {'Open Loop','Closed Loop'};
YMatrix = [(y(:,3).*sin(y(:,4))).';
           (yCL(:,3).*sin(yCL(:,4))).'];
createfigure(Nf, t, YMatrix, assi, legenda, tick, LT, LC)

Nf = 3;
assi = {'$t$ [s]','$\tilde{\omega}$ [deg/s]','Yaw rate'};
legenda = {'Open Loop','Closed Loop','Reference'};
YMatrix = [(y(:,5)*180/pi).';
           (yCL(:,5)*180/pi).';
           w_ref*180/pi];
createfigure(Nf, t, YMatrix, assi, legenda, tick, LT, LC)

Nf = 2;
assi = {'$p_x$ [m]','$p_y$ [m]','Trajectory'};
legenda = {'Open Loop','Closed Loop'};
XMatrix = [y(:,1).';
           yCL(:,1).'];
YMatrix = [y(:,2).';
           yCL(:,2).'];
createfigure(Nf, XMatrix, YMatrix, assi, legenda, tick, LT, LC)

Nf = 1;
assi = {'$t$ [s]','$\eta$ [rad]','Controller State'};
legenda = [];
YMatrix = Ce*(yCL(:,7:8).');
createfigure(Nf, t.', YMatrix, assi, legenda, tick, LT, LC)

Nf = 1;
assi = {'$t$ [s]','$\tilde{u}_{eq}$ [Nm]','Equivalent Control'};
legenda = [];
YMatrix = tau_eq;
createfigure(Nf, t.', YMatrix, assi, legenda, tick, LT, LC)

Nf = 1;
assi = {'$t$ [s]','$\tilde{\delta}_f$ [deg]','Steering Wheel Angle'};
legenda = [];
YMatrix = delta_f.*180/pi;
createfigure(Nf, t.', YMatrix, assi, legenda, tick, LT, LC)