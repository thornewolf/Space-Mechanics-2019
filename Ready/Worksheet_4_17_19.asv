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

% find FPA1
r1p = 1.7*R('Earth');
v1p = 7.1996;
a1p = 18377;
e1p = 0.46850;
p1p = a1p*(1-e1p^2);
ta1p = 46.432;
FPA1p = atan2d(r1p*e1p*sind(ta1p),p1p)

% find dv1 and alpha1
r1m = r1p; % no position change on impulse maneuver
v1m = sqrt(MU('Earth')/r1m); % circular velocity

dv1 = sqrt(v1p^2 + v1m^2 - 2*v1p*v1m*cosd(FPA1p));

alpha1 = 180 - acosd( (v1p^2 - v1m^2 - dv1^2)/(-2*v1m*dv1) );

% find FPA2

r2m = 20460 + R('Earth');
v2m = 2.8309;
a2m = a1p;
e2m = e1p;
p2m = a2m*(1-e2m^2)
ta2m = 186.43;
FPA2m = atan2d(r2m*e2m*sind(ta2m),p2m);

% find dv2 and alpha2
r2p = r2m;
v2p = sqrt(MU('Earth')/r2p); % circular velocity
FPA2p = 0;

dv2 = sqrt(v2p^2 + v2m^2 - 2*v2p*v2m*cosd(FPA2m - FPA2p))

alpha2 = 180 - acosd( (v2p^2 - v2m^2 - dv2^2)/(-2*v2m*dv2) )

% find phase angle compared to a spacecraft in the final orbit
TOF = 12375;
TA = 140;
mean_motion = sqrt(MU('Earth')/a2m^3);

syms phi
eq = TA - phi == mean_motion*TOF
phi = double(solve(eq, phi))