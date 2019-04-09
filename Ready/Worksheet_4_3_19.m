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

MU_europa = 3202.7;
a_europa = 671200;
R_europa = 1561.0;
e_europa = 0.00938;

a_clipper = 21.6*R('Jupiter');
e_clipper = 0.5731;

% 2. Europa conditions
r_europa1 = 664904.144;
v_europa1 = sqrt(2*MU('Jupiter')/r_europa1 - MU('Jupiter')/a_europa)
FPA_europa = 0; %periapsis

% 4. inf vel rel to europa
v_clipper = 17.294;
FPA_clipper = -4.5205;
v_inf = sqrt(v_europa1^2 + v_clipper^2 - 2*v_europa1*v_clipper*cosd(FPA_clipper));

% 5. fly by ang
% syms a_h_clipper
% eq = v_inf^2/2 == -MU('Jupiter')/(2*a_h_clipper)
% a_h_clipper = solve(eq, a_h_clipper)
% r_p_clipper = 
% e_h_clupper = 1 - r_p_clipper/a_h_clipper
% d = 2*asin(1/e_h_clipper)
syms FBA
dv = 1.1;
eq = dv^2 == v_inf^2 + v_inf^2 - 2*v_inf*v_inf*cosd(FBA)
FBA = max(double(solve(eq, FBA)))
