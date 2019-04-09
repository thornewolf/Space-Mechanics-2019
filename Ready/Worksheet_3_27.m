R_earth = 6378;
r_earth_sc = R_earth + 210;
r_sun_earth = 149.6e6;
v_earth = 29.785;
vt1 = 40.082;

% r1 = 0;
% a1 = 0;
% e1 = 0;
% ta1 = 0;
% 
% r2 = 0;
% a2 = 0;
% e2 = 0;
% ta2 = 0;
% 
% at = 0;
% et = 0;
% 
% syms at et
% eq(1) = r1 == at*(1-et);
% eq(2) = r2 == at*(1-et^2)/(1+et*cosd(ta2));
% [at et] = solve(eq, [at, et])

% 3.
v_infinity_sc = vt1 - v_earth;

MU_earth = 398600;
MU_sun = 132712440000;

% 4.
ah = -MU_earth/v_infinity_sc^2;
eh = 1-r_earth_sc/ah

% 5.
vc = sqrt(MU_earth/(r_earth_sc));
dv = sqrt(v_infinity_sc^2 + 2*MU_earth/r_earth_sc) - vc

% 6.
r_earth_sc = 1000 + R_earth;
vc = sqrt(MU_earth/(r_earth_sc));
dv = sqrt(v_infinity_sc^2 + 2*MU_earth/r_earth_sc) - vc