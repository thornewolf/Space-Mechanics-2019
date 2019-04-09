clc; clear;
% Constants
planets = {'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptine', 'Pluto'};
rad_list = [695990.0, 1739.2, 2439.7, 6051.9, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 25269.0, 1162.0];
mu_list = [132712440000.0, 4902.8, 22032.0, 324860.0, 398600.0, 42828.0, 126713000.0, 37941000.0, 5794500.0, 6836500.0, 981.6];
sma_list = [0.0, 384400.0, 57910000.0, 108210000.0, 149600000.0, 227920000.0, 778570000.0, 1433530000.0, 2872460000.0, 4495060000.0, 5906380000.0];
R = containers.Map(planets,rad_list);
MU = containers.Map(planets,mu_list);
r = containers.Map(planets,sma_list);
%%%%%%%

r_earth = 1.4851e8;
r_merc = 5.1818e7;
r_sc2 = r_merc;

at = 1.0015e8;
syms vt2
eq = -MU('Sun')/(2*at) == vt2^2/2 - MU('Sun')/r_sc2
vt2 = max(double(solve(eq,vt2)))

at = 1.0015e8;
et = 0.5059;
p = at*(1-et^2)
true_a2 = 30;
FPA2 = atand(r_sc2*et/p*sind(true_a2))

a_merc = 57910000;
e_merc = 0.2056;
p_merc = a_merc*(1-e_merc^2)
true_a_merc = 290;
FPA_merc = atand(r_merc*e_merc/p_merc*sind(true_a_merc))

dFPA = FPA2 - FPA_merc

syms v_merc
% eq = -MU('Sun')/(2*a_merc) == v_merc^2/2 - MU('Sun')/r_merc
v_merc = 2*sqrt(-MU('Sun')/(2*a_merc) + MU('Sun')/r_merc)
% v_merc = max(double(solve(eq,vt2)))
v_inf = sqrt(vt2^2 + v_merc^2 - 2*vt2*v_merc*cosd(dFPA))