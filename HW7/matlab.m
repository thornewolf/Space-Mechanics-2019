clc;

MU = 398600; % km^3/sec^2

rm = 9478; % km
vm = 4.961; % km/s
FPAm = 17.88; % deg

dv = 1.6; % km/s
alpha = -50; % deg

% 1. What is the initial Jason orbit semi-major axis, eccentricity, and
% true anomaly?
syms am
eq = -MU/(2*am) == vm^2/2 - MU/rm;
am = double(solve(eq, am)) % 6.6993e3 km

hm = rm*vm*cosd(FPAm);
pm = hm^2/MU;
em = sqrt(1-pm/am) % 0.5001

true_am = atan2d(rm*vm^2/MU*cosd(FPAm)*sind(FPAm),rm*vm^2/MU*cosd(FPAm)^2-1) % 160.0052


% 2. Draw the vector diagram

% 3. Determine the velocity and flight path angle immediately following the
% maneuver. Make sure to justify if dFPA is positive or negative.
syms vp FPAp
vp = sqrt(vm^2 + dv^2 - 2*vm*dv*cosd(180 - abs(alpha))) % 6.1136 km/s
dFPA = acosd(((dv^2)-(vm^2)-(vp^2))/(-2*vm*vp)) %11.5652 deg | positive because the dv is oriented towards earth (indicated by the negative alpha)

% 4. Determine the orbital characteristics following the maneuver: ap, ep,
% true_ap, dAOP
FPAp = FPAm - dFPA;
syms ap
rp = rm;
eq = -MU/(2*ap) == vp^2/2 - MU/rp;
ap = double(solve(eq, ap)) % 8.5290e+03 km

ep = sqrt((rp*vp^2/MU - 1)^2*cosd(FPAp)^2 + sind(FPAp)^2) % 0.1560

true_ap = atan2d(rp*vp^2/MU*cosd(FPAp)*sind(FPAp),rp*vp^2/MU*cosd(FPAp)^2 - 1)

true_ap = atan2d(rp*vp^2/MU*cosd(FPAp)*sind(FPAp),rp*vp^2/MU*cosd(FPAp)^2-1); % 141.4710
dAOP = -(true_ap - true_am) % 18.5342

% 5. Use GMAT to plot the two orbits. Propagate both orbits at least one
% period. Draw the vector diagram (vm, vp, alpha, dFPA, dv) and position
% vector for the maneuver on the priintout. check ap, ep, true_ap, vp.

% What would happen if you failed to perform this maneuver?