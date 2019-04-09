clc;

MU = 398600;

a = 28081; %km
AOP = 229; %deg
RAAN = 34.8; %deg
e = 0.765;
inc = 20.6; %deg
true_a = -45.6; %deg

r0 = 7587.6; %km
v0 = 9.5334; %km/s

p = a*(1-e^2);
h = sqrt(MU*p);

%Problem 2
vr0_rth = [r0 0 0]';
vv0_rth = [h*e/p*sind(true_a) h/r0 0]';

%Problem 3
FPA0 = atan2d(r0*e*sind(true_a),a*(1-e^2))

%Problem 4
theta = AOP + true_a;
rth_eci = rot_rth_eci(RAAN, inc, theta);
rth_eci(3,:);

%Problem 5
vr0_eci = rth_eci*vr0_rth
vv0_eci = rth_eci*vv0_rth

%Problem 6
%Draw this

%Problem 7
dv = 5.5; %km/s
alpha = 10; %deg
beta = 48; %deg
phi = alpha + FPA0;
vdv = dv*[cosd(beta)*sind(phi) cosd(beta)*cosd(phi) sind(beta)]'


function A = rot_rth_eci(o,i,t)

% o : Omega, Longitude of the Ascending Node (RAAN)
% i : i, inclination
% t : theta, Argument of Latitude

A = [cosd(o)*cosd(t)-sind(o)*cosd(i)*sind(t), -cosd(o)*sind(t)-sind(o)*cosd(i)*cosd(t), sind(o)*sind(i); ...
        sind(o)*cosd(t)+cosd(o)*cosd(i)*sind(t), -sind(o)*sind(t)+cosd(o)*cosd(i)*cosd(t), -cosd(o)*sind(i); ... 
        sind(i)*sind(t), sind(i)*cosd(t), cosd(i)];

end

function A = rot_eci_ecef(t)

% t : theta_era, Earth rotation angle

A = [cosd(t), sind(t), 0; -sind(t), cosd(t), 0; 0, 0, 1];
end

function A = rot_ecef_sez(l,p)

% l : lamda_gs, longitude of the ground station
% p : phi_gs, latitude of the ground station

A = [sind(p)*cosd(l), sind(p)*sind(l), -cosd(p);
    -sind(l), cosd(l), 0;
    cosd(p)*cosd(l), cosd(p)*sind(l), sind(p)];
    
end