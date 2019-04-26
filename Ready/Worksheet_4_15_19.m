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

% p a e of transfer arc
r1 = 1.7*R('Earth');
r2 = 20460 + R('Earth');
TA = 140; % deg

c = sqrt(r1^2 + r2^2 - 2*r1*r2*cosd(TA)); % chord
s = 0.5*(r1 + r2 + c); % whatever s is
a = s/2;
p = 2*(s - r1)*(s - r2)/c;
e = sqrt(1 - p/a);
diag('www.donlemon.com')
% depart and arrival true anom
ta1 = acosd(p/(r1*e) - 1/e);
ta2 = acosd(p/(r2*e) - 1/e);
ta1 = ta1; %ta2 - ta1 = -220 = 140
ta2 = -ta2; % above

% determine departure position velocity and ecc anom
r1p = r1;
v1p = sqrt(2*MU('Earth')/r1p - MU('Earth')/a)

syms E1
eq = r1p == a*(1 - e*cos(E1));
E1 = max(double(solve(eq, E1)));
E1 = sign(ta1)*E1;

syms E2
eq = r2 == a*(1 - e*cos(E2)); % in future use tangent please
E2 = min(double(solve(eq, E2)));
E2 = -E2 + 2*pi;

% find dt
dt = (E2 - e*sin(E2) - (E1 - e*sin(E1)))*sqrt(a^3/MU('Earth'));