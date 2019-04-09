clc; clear;
% Constants
planets = {'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptine', 'Pluto'};
rad_list = [695990.0, 1739.2, 2439.7, 6051.9, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 25269.0, 1162.0];
mu_list = [132712440000.0, 4902.8, 22032.0, 324860.0, 398600.0, 42828.0, 126713000.0, 37941000.0, 5794500.0, 6836500.0, 981.6];
sma_list = [0.0, 384400.0, 57910000.0, 108210000.0, 149600000.0, 227920000.0, 778570000.0, 1433530000.0, 2872460000.0, 4495060000.0, 5906380000.0];
R = containers.Map(planets,rad_list);
MU = containers.Map(planets,mu_list);
r = containers.Map(planets,sma_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Determine Earth and InSight heliocentric properties (r,v GPA) at
% departure
r_earth = r('Earth');
r_mars = r('Mars');

v_earth = sqrt(MU('Sun')/r_earth);
FPA_earth = 0;

r_sc1 = r_earth;
FPA_sc1 = 0;

r_sc2 = r_mars;
a_t = (r_sc1 + r_sc2)/2;

syms v_sc1
eq = -MU('Sun')/(2*a_t) == v_sc1^2/2 - MU('Sun')/r_sc1;
v_sc1 = max(double(solve(eq, v_sc1)));
% v_sc1 = sqrt(2*MU('Sun')/r_sc1 - MU('Sun')/a_t)

% 4. Determine v_inf_earth

v_inf_earth = v_sc1 - v_earth;

% 5. Required dv to depart earth
energy_inf_earth = v_inf_earth^2/2;
r_sc_earth = 300 + R('Earth');
% dv = sqrt(v_inf_earth^2 + 2*MU('Earth')/r_sc_earth) - sqrt(MU('Earth')/r_sc_earth);
syms v_earth_peri
eq = energy_inf_earth == v_earth_peri^2/2 - MU('Earth')/r_sc_earth;
v_earth_peri = max(double(solve(eq, v_earth_peri)));
v_earth_circ = sqrt(MU('Earth')/r_sc_earth);
dv_escape = v_earth_peri - v_earth_circ;

% 6. Mars and insight properties at arrival (r, v, FPA)
r_mars = r('Mars');
v_mars = sqrt(MU('Sun')/r_mars);
FPA_mars = 0;

r_sc2 = r_mars;
syms v_sc2
eq = -MU('Sun')/(2*a_t) == v_sc2^2/2 - MU('Sun')/r_sc2;
v_sc2 = max(double(solve(eq, v_sc2)));
FPA_sc2 = 0;

% determine v_inf_mars
v_inf_mars = abs(v_sc2 - v_mars); % Positive by convention
energy_inf_mars = v_inf_mars^2/2;
r_sc_mars = 150 + R('Mars');

syms v_mars_peri
eq = energy_inf_mars == v_mars_peri^2/2 - MU('Mars')/r_sc_mars;
v_mars_peri = max(double(solve(eq, v_mars_peri)));
v_mars_circ = sqrt(MU('Mars')/r_sc_mars);
dv_capture = v_mars_circ - v_mars_peri;
