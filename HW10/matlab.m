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

% Europa conditions
r_europa = 664904.144;
v_europa = sqrt(2*MU('Jupiter')/r_europa - MU('Jupiter')/a_europa);
FPA_europa = 0; %periapsis


% 2. Find ta_clipper1m
periapsis_europa = a_europa*(1-e_europa);
r_clipper1 = periapsis_europa;

v_clipper1m = 17.294;
FPA_clipper1m = -4.5205; %deg | from worksheet 4-3-19

ta_clipper1m = atan2d( r_clipper1*v_clipper1m^2/MU('Jupiter') * cosd(FPA_clipper1m)*sind(FPA_clipper1m),...
    r_clipper1*v_clipper1m^2/MU('Jupiter') * cosd(FPA_clipper1m) - 1);

% 3. determine periapsis_clipper_europa, alpha1
delta = 17.4; % deg | from worksheet 4-3-19

% inf vel rel to europa
v_inf = sqrt(v_europa^2 + v_clipper1m^2 - 2*v_europa*v_clipper1m*cosd(FPA_clipper1m));
a_clipper_europa1m = -MU_europa/v_inf^2;
e_clipper_europa1m = 1/sind(delta/2);
periapsis_clipper_europa1m = a_clipper_europa1m * (1-e_clipper_europa1m);

dv = 1.1; % km/s | from problem statement
nu = acosd( dv/(2*v_inf) ); % positive by inspection
offset_angle = acosd( (v_europa^2 - v_inf^2 - v_clipper1m^2)/(-2*v_inf*v_clipper1m) ); % positive by inspection
alpha = 180 - nu + offset_angle;
if alpha > 0
    alpha = -alpha; %we are told alpha is negative
end

% 4. Does the spacecraft lose energy?
% Geometrically, the spacecraft must lose energy. Observing `delta` it can
% be determined that the spacecraft must have lost energy since `delta`
% directs the dv vector opposing the initial velocity

% 5. Ahead or behind
% Observing the same diagram from problem one we can see that the v_inf
% after the gravity assist is directed more towards jupiter. This implies
% that it is an ahead pass. It is not necessary to draw a diagram in the
% europa centered view because all information is included in the jupiter
% centered view.

% 6. find r_clipper2 v_clipper2 FPA_clipper2 ta_clipper2 a_clipper2
% e_clipper2 dAOP
r_clipper2 = r_clipper1;
new_nu = 180 - alpha;
v_clipper2 = sqrt( dv^2 + v_clipper1m^2 - 2*dv*v_clipper1m*cosd(new_nu) );
FPA_clipper2 = acosd( (v_inf^2 - v_europa^2 - v_clipper2^2)/(-2*v_europa*v_clipper2) );
% dFPA = FPA_clipper2 - FPA_clipper1m not correct

FPA_clipper2 = abs(FPA_clipper2)*sign(alpha);

ta_clipper2 = atan2d(r_clipper2*v_clipper2^2/MU('Jupiter')*cosd(FPA_clipper2)*sind(FPA_clipper2),...
    r_clipper2*v_clipper2^2/MU('Jupiter')*cosd(FPA_clipper2)-1);
if sign(FPA_clipper2) > 0
    ta_clipper2 = ta_clipper2 + 180;
end

energy = v_clipper2^2/2 - MU('Jupiter')/r_clipper2;
a_clipper2 = -MU('Jupiter')/(2*energy);

e_clipper2 = sqrt( (r_clipper2*v_clipper2^2/MU('Jupiter') - 1)^2*cosd(FPA_clipper2)^2 + sind(FPA_clipper2)^2 );

dAOP = ta_clipper2 - ta_clipper1m;

