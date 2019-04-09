clc;
% Constants
planets = {'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptine', 'Pluto'};
rad_list = [695990.0, 1739.2, 2439.7, 6051.9, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 25269.0, 1162.0];
mu_list = [132712440000.0, 4902.8, 22032.0, 324860.0, 398600.0, 42828.0, 126713000.0, 37941000.0, 5794500.0, 6836500.0, 981.6];
sma_list = [0.0, 384400.0, 57910000.0, 108210000.0, 149600000.0, 227920000.0, 778570000.0, 1433530000.0, 2872460000.0, 4495060000.0, 5906380000.0];
R = containers.Map(planets,rad_list);
MU = containers.Map(planets,mu_list);
r = containers.Map(planets,sma_list);

% find vinf needed for transfer
vp = 0;
v_earth = 0;
vinf = vp - v_earth;

% Find dv needed to jump from parking orbit to hyperbola
rp = 0;

vc = sqrt(MU('Earth'));
dv = sqrt(vinf^2 - 2*MU('Earth')/rp) - vc;

% Patched conics example
rp = r('Earth');
alt = 500;
true_a2 = 140;

syms at et
et = (1 - r('Mars')/r('Earth'))/(r('Mars')/r('Earth')*cosd(true_a2) - 1);
at = r('Earth')/(1 - et);

% at earth, earth has the following
v_earth = sqrt(MU('Sun')/r('Earth'));
FPA = 0;

% s/c has the following
r_sc = r('Earth');
vp = sqrt(2*MU('Sun')/r_sc - MU('Sun')/at);
FPAp  = 0; % since as periapsis

v_infp = vp - v_earth; % This equation is only valid because FPAp is 0

% 2 body problem near mars
syms v_inf_mars vm v_mars
eq(1) = v_inf_mars == vm - v_mars;

% for circular capture
syms vc dv2 rf
eq(2) = v_inf_mars^2/2 == (vc + dv2)^2/2 - MU('Mars')/rf;
v_inf_mars = solve(eq(2), v_inf_mars);

% Patched conics example - capture
syms at
at = 1.9725e8;
vm = sqrt(2*MU('Mars')/227920000 - MU('Sun')/at);