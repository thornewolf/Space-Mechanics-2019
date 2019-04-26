clc; clear;
% Constants
planets = {'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptine', 'Pluto'};
rad_list = [695990.0, 1739.2, 2439.7, 6051.9, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 25269.0, 1162.0];
mu_list = [132712440000.0, 4902.8, 22032.0, 324860.0, 398600.0, 42828.0, 126713000.0, 37941000.0, 5794500.0, 6836500.0, 981.6];
sma_list = [0.0, 384400.0, 57910000.0, 108210000.0, 149600000.0, 227920000.0, 778570000.0, 1433530000.0, 2872460000.0, 4495060000.0, 5906380000.0];
R = containers.Map(planets,rad_list);
MU = containers.Map(planets,mu_list);
r = containers.Map(planets,sma_list);

% 2.
r('Saturn');
v_saturn = sqrt(MU('Sun')/r('Saturn'));
FPA = 0;

% 3.
FPAm = 0;
rm  = r('Saturn');
at = (r('Earth')+r('Saturn'))/2;
vp = sqrt(2*MU('Sun')/r('Saturn') - MU('Sun')/at);

% 5.
v_inf_saturn = vp - v_saturn;
v_inf_saturn = -v_inf_saturn;

% 6.
rp = R('Saturn') + 19980;
vp = sqrt(v_inf_saturn^2 + 2*MU('Saturn')/rp)
syms a
eq = -MU('Saturn')/(2*a) == vp^2/2 - MU('Saturn')/rp
a = double(solve(eq, a))

% 7.
vc = sqrt(MU('Saturn')/rp)
dv = sqrt(2*MU('Saturn')/rp - MU('Saturn')/a) - vc

% 8.
% No. It is massive