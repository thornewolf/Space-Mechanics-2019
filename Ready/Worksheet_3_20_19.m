clc;
R_earth = 6378;
MU = 398600;

r1 = 7917.2;
a1 = R_earth + 4000;
e1 = 0.3;

v1 = sqrt(2*MU/r1 - MU/a1)

syms FPA h
p = a1*(1-e1^2);
h = sqrt(MU*p);
eq = h == r1*v1*cos(FPA);
FPA = rad2deg(double(solve(eq, FPA)));
FPA1 = max(FPA); % true_a1 > 0 so FPA > 0

rt = r1;
at = 17418;
et = 0.54544;


vt  = sqrt(2*MU/r1 - MU/(at))

syms FPA h
p = at*(1-et^2);
h = sqrt(MU*p);
eq = h == rt*vt*cos(FPA);
FPA = rad2deg(double(solve(eq, FPA)));
FPAt = 0; % 0 Because definition of transfer is at periapsis

dv = sqrt(v1^2 + vt^2 - 2*v1*vt*cosd(FPAt - FPA1))

alpha = 180 - acosd(-(vt^2-v1^2 - dv^2)/(2*v1*dv))
alpha = -alpha; %alpha is negative because dv points towards the earth


