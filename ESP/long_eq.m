function f = long_eq(x)

global rho S v0 Cx0 FzW CR

kind = 1;

f = 1/2*rho*S*v0^2*Cx0+FzW.'*[mu_long(kind,x)+CR
                        mu_long(kind,x)+CR
                        CR
                        CR];
end