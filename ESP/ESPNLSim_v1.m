function dx = ESPNLSim_v1(t,x)

global m rho S g rxf rxr ryf rz  Iz theta
global Cx0 Cyb CR rk pinvH
global lambda0 v0 y0
global tmin tmax
global CL ICE Klqr Ce
global kest

rx = [rxf rxf rxr rxr];
ry = [ryf -ryf -ryf ryf];

% px = x(1);
% py = x(2);
v    = x(3);
beta = x(4);
w    = x(5);
psi  = x(6);
eta  = x(7:8);
hat_x = x(9:10);

FaB = 1/2*rho*S*v^2*[Cx0*cos(beta); Cyb*sin(beta)];

delta = SteeringWheel(t);

if CL == 0
    LL = lambda0;
    tau_eq = 0;
    track_error = zeros(2,1);
else
    LL = lambda0;
    xref = ReferenceGenerator(delta);
    track_error = [v*sin(beta); w]-xref;
%     veq = TrajectoryTrackingControl(xref,zeros(2,1));
    tau_eq = -Klqr(2:3)*([v*sin(beta);w]-xref)-Klqr(1)*Ce*eta;
end


tauBi = zeros(1,4);
muBi  = zeros(2,4);
for i = 1:4
    deltai = delta(i);
    
    ryi = ry(i);
    rxi = rx(i);
    VWi = [cos(deltai) sin(deltai);
        -sin(deltai) cos(deltai)]*(w*[-ryi; rxi]+v*[cos(beta); sin(beta)]);
    Vxi = VWi(1);
    Vyi = VWi(2);
    betai = atan(Vyi/Vxi);
       
    kind = rk;
    if i == 3 || i == 4
        if (ICE == 1) && (t > max(tmin)) && (t < max(tmax))
            kind = 6;
        end
    end
    
    muL = mu_long(kind,LL(i));
    muS = -mu_long(kind,betai/(pi/2));
    muB = [cos(deltai) -sin(deltai)
        sin(deltai) cos(deltai)]*[muL+CR;
        muS];
    
    muBi(:,i) = muB;
    tauBi(i) = [(sin(deltai)*rxi - cos(deltai)*ryi)   (cos(deltai)*rxi+sin(deltai)*ryi)]*muB;
end

% b = [m*g; 0; 0];
% arm_x = zeros(1,4);
% arm_y = zeros(1,4);
% for i = 1:4
%     arm_x(i) = ry(i)+rz*muBi(2,i);
%     arm_y(i) = -rx(i)-rz*muBi(1,i);
% end
% H = [1 1 1 1;
%     arm_x;
%     arm_y];
% FzW=(H.'/(H*H.'))*b;
FzW = pinvH*(rz/2*rho*S*v^2*[0; Cyb*sin(beta); -Cx0*cos(beta)]+[m*g*cos(theta); -rz*m*v*w*cos(beta); -m*g*sin(theta)-rz*m*v*w*sin(beta)]);

FB = FaB;
tauB = 0;
for i = 1:4
    FB = FB + muBi(:,i)*FzW(i) + 100*[v0-v; 0];
    tauB = tauB + tauBi(i)*FzW(i)  + tau_eq;
end

% plant
dp = [v*cos(psi+beta);
    v*sin(psi+beta)];
dV = [m^(-1) 0         0;
    0     (m*v)^(-1) 0;
    0      0         Iz^(-1)]*([cos(beta) sin(beta) 0;
    -sin(beta) cos(beta) 0;
    0         0      1]*[FB; tauB]-[0; m*v; 0]*w);
dpsi = w;

% observer
y = Sensors(x,t);
deltaf = (delta(1)+delta(2))/2;
utilde = [deltaf; tau_eq];
ytilde = (y-y0);
dot_hat_x = kest.A*hat_x + kest.B*[utilde; ytilde];
%% overall system = plant + IM + observer

dx = [dp;
    dV.*[0;1;1];
    dpsi;
    track_error;
    dot_hat_x];