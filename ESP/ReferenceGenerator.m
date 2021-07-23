function xref = ReferenceGenerator(delta)
%disp('LINEAR EQUIIBRIUM POINT')

global v0 

delta_f = (delta(1)+delta(2))/2;

g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
Cx0 = -0.3; % [-] drag coefficient
Cyb = -0.3; % [-] drag coefficient

m = 1200; % [kg] vehicle mass
rxf = 1.4;
rxr = -1.6;
ryf = 0.90;
Iz = 1200*(rxr^2+ryf^2)^2;
rx = [rxf rxf rxr rxr];
ry = [ryf -ryf -ryf ryf];
rz =  0.70;

H = [ones(1,4);
    ry;
    -rx];


beta0 = 0;
w0 = 0;
theta = 0/180*pi;

pinvH = pinv(H);
FzW = pinvH*(rz/2*rho*S*v0^2*[0; Cyb*sin(beta0); -Cx0*cos(beta0)]+[m*g*cos(theta); -rz*m*v0*w0*cos(beta0); -rz*m*g*sin(theta)-rz*m*v0*w0*sin(beta0)]);
Cf = (FzW(1)+FzW(2))*Partial_mu_long(1,0);
Cr = (FzW(3)+FzW(4))*Partial_mu_long(1,0);

A11 = -(Cf+Cr)/(m*v0);
A12 = -v0-(Cf*rxf+Cr*rxr)/(m*v0);
A21 = -(Cf*rxf+Cr*rxr)/(Iz*v0);
A22 = -(Cf*rxf^2+Cr*rxr^2)/(Iz*v0);
A = [A11 A12
     A21 A22];
B1 = Cf*[1/m;
      rxf/Iz];
xref = -inv(A)*B1*delta_f;
% if abs(xeq(2)) > 0.85*1.3*g/v0
%     xeq(2) = 0.85*1.3*g/v0 * sign(xeq(2));
% end
% if abs(xeq(1)) > atan(0.02*1.3*g)
%     xeq(1) = atan(0.02*1.3*g) * sign(xeq(1));
% end
% disp(['beta = ',num2str(xeq(1)*180/pi),' deg'])
% disp(['omega = ',num2str(xeq(2)*180/pi),' deg/s'])
end