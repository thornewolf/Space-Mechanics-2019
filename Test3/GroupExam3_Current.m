%% Group Exam 3
clear; clc;
% AE313 Exam #3  - Due: 4/10/19
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
vr1_RTHBefore = [rInit; 0; 0];
vr1_MCIBefore = A*vr1_RTHBefore;
vv1_EPHBefore = sqrt(mu('Mars')/pInit)*[ -sind(taInit); eInit + cosd(taInit); 0 ];
vv1_RTHBefore = B*vv1_EPHBefore;
vv1_MCIBefore = A*vv1_RTHBefore;
v1 = norm(vv1_MCIBefore);
fpa1 = atan2d(rInit*eInit*sind(taInit), pInit);
%%
fprintf('Problem 1:\n')
fprintf('Initial Radius in Mars Centered Inertial (MCI): %gx + %gy + %gz km \n',vr1_MCIBefore)
fprintf('Initial Velocity in Mars Centered Inertial (MCI): %gx + %gy+ %gz km \n\n',vv1_MCIBefore)

%% Problem 2

% The first burn is an in-plane burn with |deltaV1| = 0.14 km/s and alpha1
% = 45 degrees. Draw the vector diagram.

deltaV1 = 0.14; %km/s --
alpha1 = 45; %degrees --

%% Problem 3

% Find the velocity (magnitude) and flight path angle immediately after the
% first maneuver.
r1After = rInit;
v1After = sqrt(deltaV1^2 + v1^2 - 2*deltaV1*v1*cosd(180-alpha1));
dfpa = asind(sind(180-alpha1)/v1After * deltaV1);
fpaAfter = fpa1 + dfpa;

%%
fprintf('Problem 3:\n')
fprintf('Velocity magnitude: %g km/s\n',v1After)
fprintf('Flight Path Angle: %g degrees \n\n', fpaAfter)

%% Problem 4

% Wait until MAVEN reaches TA = 240 degrees in the new orbit before
% implementing the final maneuver. What is the new position and velocity
% (r2, v2) in the Mars equatorial inertial coordinate system (x - y - z)?

ta1After = atan2d(r1After*v1After^2/mu('Mars')*sind(fpaAfter)*cosd(fpaAfter), r1After*v1After^2/mu('Mars')*cosd(fpaAfter)^2-1);
ta2Before = 240; % (deg) final true anomaly
syms a2Before
a2Before = double(solve(-mu('Mars')/(2*a2Before) == v1After^2/2 - mu('Mars')/r1After, a2Before));
e2Before = double(sqrt((r1After*v1After^2/mu('Mars')-1)^2*cosd(fpaAfter)^2+sind(fpaAfter)^2));
p2Before = (a2Before)*(1 - e2Before^2);
h2Before = sqrt(mu('Mars')*p2Before);
i2Before = incInit;
raan2Before = raanInit;
r2Before = p2Before/(1+e2Before*cosd(ta2Before));
rHat2Before =  vr1_MCIBefore/norm(vr1_MCIBefore);
hHat2Before = cross(vr1_MCIBefore, vv1_MCIBefore) / norm(cross(vr1_MCIBefore, vv1_MCIBefore));
thetaHat2Before = cross(hHat2Before, rHat2Before);
aol2Before = asind(rHat2Before(3)/sind(i2Before));
aol2BeforeCheck = acosd(thetaHat2Before(3)/sind(i2Before)); % Choosing pos. value
aol2Before = aol2BeforeCheck;
deltaAOP1 = taInit - ta1After; % (deg) change in argument of periapsis of transfer
aop2After = aopInit + deltaAOP1; % (deg) argument of periapsis of transfer
aol2After = ta2Before + aop2After; % (deg) argument of latittude of transfer
A2 = RTHtoMCI(raanInit, aol2After, incInit);
B2 = EPHtoRTH(ta2Before);
vv2_RTHBefore = [h2Before*e2Before/p2Before*sind(ta2Before) h2Before/r2Before 0]';
vr2_RTHBefore = [r2Before; 0; 0];
vr2_MCIBefore = double(A2*vr2_RTHBefore);
vv2_MCIBefore = double(A2*vv2_RTHBefore);

%%
fprintf('Problem 4:\n')
fprintf('New Radius Vector: %gx + %gy + %gz km \n',vr2_MCIBefore)
fprintf('New Velocity Vector: %gx + %gy + %gz km/s \n\n',vv2_MCIBefore)

%% Problem 5

% The next maneuver has |deltaV2| = 0.22 km/s, alpha2 = -20 degrees, and
% beta2 = -50 degrees. Find deltaV2 in the Mars equitorial inertial
% coordinate system
syms a2After
deltaV2 = 0.22; % --
alpha2 = -20; % --
beta2 = -50; % --
h2Before = norm(cross(vr2_MCIBefore, vv2_MCIBefore));
fpa2Before = atan2d(r2Before*e2Before*sind(ta2Before), p2Before);
phi = alpha2 + fpa2Before;
vdeltaV2_RTH = deltaV2*[cosd(beta2)*sind(phi); cosd(beta2)*cosd(phi); sind(beta2)];
vdeltaV2_MCI = A2*vdeltaV2_RTH;

%%
fprintf('Problem 5:\n')
fprintf('deltaV Vector after Manuever 2 in MCI: %gx + %gy + %gz km/s \n\n', vdeltaV2_MCI)

%% Problem 6

% Determine the new velocity in the Mars equatorial inertial coordinate
% system.
vv2_MCIAfter = vv2_MCIBefore + vdeltaV2_MCI;


%%
fprintf('Problem 6:\n')
fprintf('New Velocity in MCI: %gx + %gy + %gz km/s \n\n',vv2_MCIAfter);

%% Problem 7

% How have the orbital elements changed from the initial orbit (values
% given in the problem statement)? What do these changes mean for the new
% orbit?
syms aFinal eFinal taFinal

rHatFinal = vr2_MCIBefore / norm(vr2_MCIBefore);
hHatFinal = cross(vr2_MCIBefore,vv2_MCIAfter)/norm(cross(vr2_MCIBefore,vv2_MCIAfter));
thetaHatFinal = cross(hHatFinal,rHatFinal);

hFinal = norm(cross(vr2_MCIBefore,vv2_MCIAfter));
aFinal = double(solve(-mu('Mars')/(2*aFinal) == norm(vv2_MCIAfter)^2/2 - mu('Mars')/norm(vr2_MCIBefore), aFinal));
eFinal = max(double(solve(-mu('Mars')^2/(2*hFinal^2) * (1-eFinal^2) == -mu('Mars')/(2*aFinal), eFinal)));
incFinal = acosd(hHatFinal(3)); % Choose positive value because 0 < i < 180
raanFinal = acosd(-(hHatFinal(2)/sind(incFinal))); 
raanFinalCheck = asind(hHatFinal(1)/sind(incFinal));
raanFinal = raanFinalCheck; % Choose the negative value
aolFinal = asind(rHatFinal(3)/sind(incFinal)); 
aolFinalCheck = acosd(thetaHatFinal(3)/sind(incFinal));
aolFinal = aolFinalCheck; % Choose the positive value
pFinal = aFinal * (1 - eFinal^2);
taFinal = -acosd(((pFinal/norm(vr2_MCIBefore)) - 1)/eFinal); % Choose negative because descending (alpha negative)
aopFinal = aolFinal - taFinal;

deltaSMA = aFinal - aInit;
deltaECC = eFinal - eInit;
deltaRAAN = 360 + (raanFinal - raanInit);
deltaAOP = aopFinal - aopInit;
deltaINC = incFinal - incInit;
deltaAOL = aolFinal - aolInit;

%%
fprintf('Problem 7:\n')
fprintf('Change in a: %f km \n', deltaSMA)
fprintf('Change in e: %f \n', deltaECC)
fprintf('Change in RAAN: %f degrees \n', deltaRAAN)
fprintf('Change in AOP: %f degrees \n', deltaAOP)
fprintf('Change in Inclination: %f degrees \n', deltaINC)
%%
% Increasing the semi-major axis scales the orbit away
% from Mars such that the spacecraft is near apoapsis for a longer period
% of time, allowing for more communication. The eccentricity is about the
% same indicating a minimal change in the elognation of the orbit. The
% change in the Right Ascention of the Ascending Node is very small. This
% means that the orbit did not rotate significantly, implying the Ascending
% Node and the Descending Node are in similar places in space. The Argument
% of Periapsis rotates the orbit's closest approach about Mars. The change
% in Inclination shows that the orbit is pulled further from the poles and
% is more of an equatorial orbit.


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

%% End