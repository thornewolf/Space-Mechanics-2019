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

et = 0.21378;
at = 190276560;

v_sc = 21.612;
v_mars = 24.130;
FPA_sc = 4.74;
v_inf = sqrt(v_sc^2 + v_mars^2 - 2*v_sc*v_mars*cosd(FPA_sc))

theta = asind(v_sc/v_inf*sind(FPA_sc));
d = theta*2

dv = v_inf*sqrt(2 - 2*cosd(d))
abs_a = 180 - (180-2*FPA_sc)/2
a = -abs_a

%last
syms a_h
eq = v_inf_merc^2/2 == -MU('Mercury')/(2*a_h);
a_h = solve(eq, a_h);
a_h = double(a_h);

r_p = 200 + R('Mercury');

e_h = 1 - r_p/a_h;

syms d
eq = sind(d/2) == 1/e_h
d = min(double(solve(eq, d)))