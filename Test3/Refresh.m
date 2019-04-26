%% Group Exam 3
clear; clc;
% AE313 - Due: 4/10/19
% Grace Day, JiYoung Hwang, Aaron Scott, Thorne Wolfenbarger

%% Preamble
bodies = { 'Earth', 'Mars', 'Sun' };
gravParams = [ 398600, 42828, 132712440000 ];
radii = [ 6378, 3397, 695990 ];
sma = [149600000, 227920000, 0];
mu = containers.Map(bodies, gravParams); % (km^3/s^2) Gravitational parameter
R = containers.Map(bodies, radii); % (km) Mean equitorial radius
a = containers.Map(bodies, sma); % (km) Semi-major axis

%% Given
aInit = 5744; % km --
eInit = 0.3842; % --
raanInit = 321.7; % degrees --
aopInit = 303.2; % degrees --
incInit = 74.77; % degrees --
aolInit = 203.9; % degrees --

%% Problem 1
% Find the current position and velocity (r1, v1) vectors in Mars
% equatorial inertial system (x - y - z), as well as the initial flight
% path angle.
taInit = aolInit - aopInit; % Initial true anomaly value
A = RTHtoMCI(raanInit, aolInit, incInit);
B = EPHtoRTH(taInit);
pInit = aInit*(1-eInit^2);
rInit = pInit/(1+eInit*cosd(taInit));
vr1_RTH = [rInit; 0; 0];
vr1_MCI = A*vr1_RTH;
vv1_EPH = sqrt(mu('Mars')/pInit)*[ -sind(taInit); eInit + cosd(taInit); 0 ];
vv1_RTH = B*vv1_EPH;
vv1_MCI = A*vv1_RTH;
v1 = norm(vv1_MCI);
fpa1 = atan2d(rInit*eInit*sind(taInit), pInit);

%% Problem 2
% The first burn is an in-plane burn with |deltaV1| = 0.14 km/s and alpha1
% = 45 degrees. Draw the vector diagram.

deltaV1 = 0.14; %km/s --
alpha1 = 45; %degrees --



%% Function Definitions
function outputMatrix = RTHtoMCI(raan, theta, i)
    outputMatrix = ...
        [ cosd(raan)*cosd(theta) - sind(raan)*cosd(i)*sind(theta) , ...
         -cosd(raan)*sind(theta) - sind(raan)*cosd(i)*cosd(theta) , ...
          sind(raan)*sind(i) ; ...
          sind(raan)*cosd(theta) + cosd(raan)*cosd(i)*sind(theta) , ...
         -sind(raan)*sind(theta) + cosd(raan)*cosd(i)*cosd(theta) , ...
         -cosd(raan)*sind(i) ; ...
          sind(i)*sind(theta) , ...
          sind(i)*cosd(theta) , ...
          cosd(i) ];
end

function outputMatrix = EPHtoRTH(ta)
    outputMatrix = ...
        [ cosd(ta), sind(ta), 0; ...
         -sind(ta), cosd(ta), 0; ...
                 0,        0, 1; ...
        ];
end