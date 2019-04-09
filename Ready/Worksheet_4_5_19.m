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

r_merc = 5.1818e7;
v_merc = 53.203;
FPA_merc = -10.232;

e_t = 0.5059;
a_t = 1.0015e8;

v_2 = 61.521;
FPA_2 = 9.9756;

% 3.
p_t = a_t*(1-e_t^2)
FPA_t = 9.9756;
dFPA = FPA_merc - FPA_t

v_inf_merc = sqrt(v_merc^2 + v_2^2 - 2*v_merc*v_2*cosd(dFPA))

syms a_h
eq = v_inf_merc^2/2 == -MU('Mercury')/(2*a_h);
a_h = solve(eq, a_h);
a_h = double(a_h);

r_p = 200 + R('Mercury');

e_h = 1 - r_p/a_h;

syms d
eq = sind(d/2) == 1/e_h
d = min(double(solve(eq, d)))