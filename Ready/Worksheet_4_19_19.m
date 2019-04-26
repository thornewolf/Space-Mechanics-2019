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

r1 = R('Earth')*1.7;
r2 = R('Earth') + 20460;
c = 35828.32;
TA = 140;

emin = (r2 - r1)/c;

xi = atan2d( (r1/r2 - cosd(TA)),sind(TA) );

% True Anomaly Range
ta_upper = atand(-tand(TA/2));
if ta_upper < 0
    ta_upper = ta_upper + 180;
end
ta_range = [-xi  ta_upper];

% find et and at
ta1p = 55; % deg
et  = emin / sind(ta1p + xi);

pt = r1 * (1 + et*cosd(ta1p));

at = pt / (1 - et^2);

% conditions on transfer arc at departure
r1p = r1;
v1p = sqrt(2*MU('Earth')/r1p - MU('Earth')/at);

FPA1p = atan2d(r1p * et * sind(ta1p), pt);

% find dv and alpha1p
vc = sqrt(MU('Earth') / r1);
v1m = vc;
dv1 = sqrt( v1m^2 + v1p^2 - 2*v1m*v1p*cosd(FPA1p) );

nu = acosd( (v1p^2 - v1m^2 - dv1^2)/(-2*v1m*dv1) );

alpha = 180 - nu;

% find conditions at arrival (v2m r2m FPA2m)
r2m = r2;
v2m = sqrt(2*MU('Earth')/r2m - MU('Earth')/at);

ta2m = ta1p + TA;
FPA2p = atan2d(r2m * et * sind(ta2m), pt)

vc = sqrt(MU('Earth')/r2);
v2p = vc;

dv2 = sqrt( v2m^2 + v2p^2 - 2*v2m*v2p*cosd(FPA2p) );

nu = acosd( (v2p^2 - v2m^2 - dv2^2)/(-2*v2m*dv2) );

alpha = 180 - nu;

dvt = dv1 + dv2;