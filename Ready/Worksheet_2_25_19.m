clc;
MU = 398600
r0 = 6761
v0 = 8.124
gamma0 = 12.78

% Question 1
syms v r a
energy = v^2/2 - MU/r == -MU/(2*a);
energy_s = subs(energy, [r v] , [r0 v0]);
a0 = double(solve(energy_s,a))

% Question 2
syms e gamma
ecc_eq = e == ((r*v^2/MU-1)^2*cosd(gamma)^2+sind(gamma)^2)^0.5;
ecc_eq_s = subs(ecc_eq, [r v gamma] , [r0 v0 gamma0]);
double(solve(ecc_eq_s,e))

% Question 3
syms true_a
true_a = atan2d((r*v^2/MU*cosd(gamma)*sind(gamma)),(r*v^2/MU*cosd(gamma)^2 -1));
true_a0 = double(subs(true_a, [r v gamma] , [r0 v0 gamma0]))

% Question 4
% Drawing

% Question 5
% Same as original position

% Question 6
v0 = 8.124;
dv = 0.9;
v1 = v0 + dv

% Question 7
% Same FPA as original
true_a = atan2d((r*v^2/MU*cosd(gamma)*sind(gamma)),(r*v^2/MU*cosd(gamma)^2 -1));
true_a1 = double(subs(true_a, [r v gamma] , [6761 v1 12.78]))

% Question 8
energy = v^2/2 - MU/r == -MU/(2*a);
energy_s = subs(energy, [r v] , [6761 v1]);
a1 = double(solve(energy_s,a))

ecc_eq = e == ((r*v^2/MU-1)^2*cosd(gamma)^2+sind(gamma)^2)^0.5;
ecc_eq_s = subs(ecc_eq, [r v gamma] , [6761 v1 12.78]);
ecc1 = double(solve(ecc_eq_s,e))