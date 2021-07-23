function veq = TrajectoryTrackingControl(xref,dotxref)

global Klqr A Iz
K = -Klqr;
B2eq = [0; 1/Iz];
Cstar = [0 1];
veq = -inv(Cstar*inv(A+B2eq*K)*B2eq)*Cstar*(xref-inv(A+B2eq*K)*dotxref);
end