function f = SteeringEq_v2(x)

global rk
global V delta a b N

g = 9.81; % [m/s^2] gravity acceleration
rho = 1.225; % [kg/m^3] air density
S = 2.4; % [m^2] cross surface
CD = 2.2; % [-] drag coefficient

M = 1200; % [kg] vehicle mass

R = x(1);
beta = x(2);
lambdaL1 = x(3);

Fb = [sin(beta) cos(beta);
      cos(beta) -sin(beta)]*[M*V^2/R; 1/2*rho*S*V^2*CD];
Fxb = Fb(1);
Fyb = Fb(2);

Fw = [cos(delta) -sin(delta) 0;
      sin(delta) cos(delta) +1;
      a*sin(delta) a*cos(delta) -b]\[Fxb; Fyb; 0];
  
Fw1x = Fw(1);
Fw1y = Fw(2);
Fw2y = Fw(3);

muL1 = Fw1x/N;
muS1 = Fw1y/N;
muS2 = Fw2y/N;

Vw1b = V*[cos(beta); sin(beta)] + [0;1]*V/R*a;
Vw1 = [cos(delta) sin(delta);
      -sin(delta) cos(delta)]*Vw1b;
w1r = Vw1(1)/(1+lambdaL1);
lambdaS1 = Vw1(2)/w1r;
lambdaT1 = sqrt(lambdaL1^2+lambdaS1^2);

Vw2b = V*[cos(beta); sin(beta)] - [0;1]*V/R*b;
Vw2 = Vw2b;
w2r = Vw2(1);
lambdaS2 = Vw2(2)/w2r;
lambdaT2 = sqrt(lambdaS2^2);

f = [muL1 + lambdaL1/lambdaT1*mu_long(rk,lambdaT1);
     muS1 + lambdaS1/lambdaT1*mu_long(rk,lambdaT1);
     muS2 + lambdaS2/lambdaT2* mu_long(rk,lambdaT2)];