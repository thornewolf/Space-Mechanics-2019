clc; clear all; close all
% 2. Find vdv in intertial unit vectors
MU = 398600; %kg
a = 28081; %km
AOP = 229; %deg
RAAN = 35.8; %deg
e = 0.765; %deg
inc = 20.6; %deg
true_a = -45.6; %deg

vr1_eci = [-5978.4 -4668.0 -158.31]';
vv1_eci = [7.7465 -4.6198 -3.0877]';

rth_eci = [-0.788, 0.5820, 0.2008;...
           -0.6153, -0.7334, -0.2889;...
           -0.02087, -0.3512, 0.9361];
       
dv = 5.5; %km/s
alpha = 10; %deg
beta = 47; %deh
vdv_rth = [-0.62533 3.6985 4.0224]';

% 2.
vdv_eci = rth_eci*vdv_rth;

% 4.
vv2_eci = vv1_eci + vdv_eci;

% 5.
FPA2 = asind(dot(vr1_eci,vv2_eci)/(norm(vr1_eci)*norm(vv2_eci)))

% 6.
h = cross(vr1_eci,vv2_eci);
h_hat = h/norm(h);
inc2 = acosd(h_hat(3)); %only defined positive
dinc = inc2 - inc;

% 7.
r2 = norm(vr1_eci)
v2 = norm(vv2_eci)
true_a2 = atan2d(r2*v2^2/MU*cosd(FPA2)*sind(FPA2),r2*v2^2/MU*cosd(FPA2)^2-1)


AOP2 = true_a2 - true_a %Wrong. Only valid for 2d maneuvers
dAOP = AOP2 - AOP






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