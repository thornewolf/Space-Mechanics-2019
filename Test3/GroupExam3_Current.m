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

%% Problem 3
% Find the velocity (magnitude) and flight path angle immediately after the
% first maneuver.

cInc = incInit; % It was a planar maneuver

% Can use same rotation matrix for same r-value. Sanity check
r1After = rInit;
vr1After_RTH = [r1After; 0; 0];
vr1After_MCI = A*vr1After_RTH;

% Get dv as a vector in eci
dv1 = deltaV1;
beta1 = 0;
vdv1_VCN = dv1*[cosd(beta1)*cosd(alpha1) cosd(beta1)*sind(alpha1) sind(beta1)]';
phi1 = fpa1 + alpha1;
vdv1_RTH = dv1*[cosd(beta1)*sind(phi1) cosd(beta1)*cosd(phi1) sind(beta1)]';

% Need to update the rotation matrix for the new orbit properties
% Get rth vectors
rhat = vr1After_MCI/norm(vr1After_MCI);
h = cross(vr1After_MCI,vv1_MCI);
hhat = h/norm(h);
that = cross(hhat,rhat);



% Calculate new RAAN value
syms cRAAN
eq(1) = hhat(1) == sind(cRAAN)*sind(cInc);
eq(2) = hhat(2) == -cosd(cRAAN)*sind(cInc);
cRAAN1 = double(solve(eq(1), cRAAN));
cRAAN2 = double(solve(eq(2), cRAAN));
% Checks the matching value
if abs(cRAAN1(1) - cRAAN2(1)) < 0.01 | abs(cRAAN1(1) - cRAAN2(2)) < 0.01
    cRAAN = cRAAN1(1);
else
    cRAAN = cRAAN1(2);
end
if cRAAN < 360
    cRAAN = cRAAN + 360;
end
syms cAOL
eq(1) = rhat(3) == sind(cInc)*sind(cAOL);
eq(2) = that(3) == sind(cInc)*cosd(cAOL);
cAOL1 = double(solve(eq(1), cAOL));
cAOL2 = double(solve(eq(2), cAOL));
if abs(cAOL1(1) - cAOL2(1)) < 0.01 | abs(cAOL1(1) - cAOL2(2)) < 0.01
    cAOL = cAOL1(1)
else
    cAOL = cAOL1(2)
end

A = RTHtoMCI(cRAAN, cAOL, cInc);

v1After = sqrt(deltaV1^2 + v1^2 - 2*deltaV1*v1*cosd(180-alpha1));
vv1After_RTH =[v1After; 0; 0];
vv1After_MCI = A*vv1After_RTH; %Incorrect but magnitude is correct
vdv1_MCI = A*vdv1_RTH; %Incorrect
altvv1After = vdv1_MCI + vv1_MCI;

dfpa = asind(sind(180-alpha1)/v1After) * deltaV1;
fpaAfter = fpa1 + dfpa;

true_a_after = 0;

%% Problem 4
% Wait until MAVEN reaches TA = 240 degrees in the new orbit before
% implementing the final maneuver. What is the new position and velocity
% (r2, v2) in the Mars equatorial inertial coordinate system (x - y - z)?
syms a2Before p 
a2Before = solve(-mu('Mars')/(2*a2Before) == v1After^2/2 - mu('Mars')/r1After, a2Before);
p2Before = solve(sqrt(mu('Mars')*p) == r1After*v1After*cosd(fpaAfter),p);
e2Before = sqrt(1 - p2Before/a2Before);
ta2Before = 240; % --
i2Before = incInit;
raan2Before = raanInit;
r2Before = p2Before/(1+e2Before*cosd(ta2Before));
v2Before = 0;
rHat2Before =  vr1_MCI/norm(vr1_MCI);
hHat2Before = cross(vr1_MCI, vv1_MCI) / norm(cross(vr1_MCI, vv1_MCI));
thetaHat2Before = cross(hHat2Before, rHat2Before);
aol2Before = asind(rHat2Before(3)/sind(i2Before))
aol2BeforeCheck = acosd(thetaHat2Before(3)/sind(i2Before)) % Choosing pos. value
aol2Before = aol2BeforeCheck;
B2 = EPHtoRTH(ta2Before);
vr2_RTHBefore = [r2Before; 0; 0];
vr2_MCIBefore = double(A*vr2_RTHBefore);
vv2_EPHBefore = sqrt(mu('Mars')/pInit)*[ -sind(ta2Before); e2Before + cosd(ta2Before); 0 ];
vv2_RTHBefore = B2*vv2_EPHBefore;
vv2_MCIBefore = double(A*vv2_RTHBefore);

% raan2Before = asind(hHat2Before(1)/sind(i2Before));
% raan2Before = [raan2Before 180-raan2Before]
% raan2BeforeCheck = acosd(-hHat2Before(2)/sind(i2Before)); % Consistent with this value
% raan2BeforeCheck = [raan2BeforeCheck -raan2BeforeCheck]
% aol2Before = 0;

%% Problem 5
% The next maneuver has |deltaV2| = 0.22 km/s, alpha2 = -20 degrees, and
% beta2 = -50 degrees. Find deltaV2 in the Mars equitorial inertial
% coordinate system

deltaV2 = 0.22; % --
alpha2 = -20; % --
beta2 = -50; % --
vdeltaV2_VCN = deltaV2*[cosd(beta2)*cosd(alpha2); cosd(beta2)*sind(alpha2); sind(beta2)];
vdeltaV2_MCI = A*vdeltaV2_VCN

%% Problem 6
% Determine the new velocity in the Mars equatorial inertial coordinate
% system.


%% Problem 7
% How have the orbital elements changed from the initial orbit (values
% given in the problem statement)? What do these changes mean for the new
% orbit?

%% Values for GMAT Check
vdeltaV1_VCN = deltaV1*[cosd(alpha1); sind(alpha1); 0]; % km/s - Mars Centered Inertial deltaV
vdeltaV2_VCN = deltaV2*[cosd(beta2)*cosd(alpha2); cosd(beta2)*sind(alpha2); sind(beta2)]; % km/s - Mars Centered Inertial deltaV

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