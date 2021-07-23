clc
close all
clear all

m = 1600; %[kg] mass
g = 9.81; % [m/s^2] gravity acceleration

r = 1;
a = 1.4;
b = 1.6;
d = 0.5;
L = a+b;

H = [-2*r -2*r 0 0
    L 0 L 0
    1 1 1 1];
M = (H'/(H*H'));
R = 1000;

V = linspace(0,37*3.6);
ar = (V.^2)./R;

FO = m/(2*r)*(ar*d+g*r);
FI = m/(2*r)*(-ar*d+g*r);

FR = m*g*a/L;
FF = m*g*b/L;

H2 = [ 1 1 0 0
    1 0 1 0
    1 1 1 1];



for i = 1:length(ar)
%     Y = [ar(i)*d-m*g*r;
%         b*m*g;
%         m*g];
%     
%     X = M*Y;
%     FOF(i) = X(1);
%     FOR(i) = X(2);
%     FIF(i) = X(3);
%     FIR(i) = X(4);
        
    Y2 = [FI(i);
         FR;
         m*g];
    
    X2 = (H2'/(H2*H2'))*Y2;
    FIR(i) = X2(1);
    FIF(i) = X2(2);
    FOR(i) = X2(3);
    FOF(i) = X2(4);
    
end

figure
subplot(2,2,1)
plot(ar/g,FOF./(m.*g))
grid on
xlabel('$a_R/g\,\, \left[-\right]$','Interpret','latex')
title('$F_z^{OF}/ mg$','Interpret','latex')
xlim([0 max(ar/g)])

subplot(2,2,3)
plot(ar/g,FOR./(m.*g))
grid on
xlabel('$a_R/g\,\, \left[-\right]$','Interpret','latex')
title('$F_z^{OR}/ mg$','Interpret','latex')
xlim([0 max(ar/g)])

subplot(2,2,2)
plot(ar/g,FIF./(m.*g))
grid on
xlabel('$a_R/g\,\, \left[-\right]$','Interpret','latex')
title('$F_z^{IF}/ mg$','Interpret','latex')
xlim([0 max(ar/g)])

subplot(2,2,4)
plot(ar/g,FIR./(m.*g))
grid on
xlabel('$a_R/g\,\, \left[-\right]$','Interpret','latex')
title('$F_z^{IR}/ mg$','Interpret','latex')
xlim([0 max(ar/g)])
