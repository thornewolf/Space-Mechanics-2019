clc;

r2_m = 14179; % km
r2_p = r2_m
v2_m = 4.4434; % km/s
FPA_m = -19.2; % deg
dv = 0.5; % km/s
a2_m = 10927; % km
e2_m = 0.4326;
alpha = -15; % deg
MU = 1;

n = 180 - abs(alpha);

syms v2_p
eq = v2_p^2 == dv^2 + v2_m^2 - 2*dv*v2_m*cosd(n);
v2_p = double(solve(eq,v2_p));
v2_p = max(v2_p); % Answer for 4

syms dgamma
eq = dv^2 == v2_p^2 + v2_m^2 - 2*v2_p*v2_m*cosd(dgamma);
dgamma = min(double(solve(eq,dgamma)))

syms FPA_p
eq = dgamma + FPA_m == FPA_p;
FPA_p = double(solve(eq, FPA_p))

syms true_a_p
eq = true_a_p == atan2d(r2_p*v2_p^2/MU*cosd(dgamma)*sind(dgamma),r2_p*v2_p^2/MU*cosd(dgamma)^2-1)
true_a_p = double(solve(eq, true_a_p)) % Wrong