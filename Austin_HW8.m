% Austin Kendall
% Homework 8
% 29 March 2019
clc; clear all

%% Constants

mu = 398600;
a1 = 28081;
w1 = 229;
Omega1 = 34.8;
e1 = 0.765;
i1 = 20.6;
thetastar1 = -45.6;
p1 = a1*(1-e1^2);
theta1 = thetastar1 + w1;

rminus_vec = [-5978.4;-4668;-158.31];
rminus = norm(rminus_vec);
vminus_vec = [7.7465;-4.6198;-3.0877];
vminus = norm(vminus_vec);

gammaminus = atan2d(rminus*e1*sind(thetastar1),p1);

dv1 = 4.5;
alpha1 = -30;
beta1 = 40;

%% 1
phi1 = alpha1 + gammaminus;
dv1_rth = dv1*[cosd(beta1)*sind(phi1);cosd(beta1)*cosd(phi1);sind(beta1)];

A_rth_xyz = [cosd(Omega1)*cosd(theta1)-sind(Omega1)*cosd(i1)*sind(theta1), ...
    -cosd(Omega1)*sind(theta1)-sind(Omega1)*cosd(i1)*cosd(theta1), ... 
        sind(Omega1)*sind(i1); sind(Omega1)*cosd(theta1)+cosd(Omega1)*cosd(i1)*sind(theta1), ...
            -sind(Omega1)*sind(theta1)+cosd(Omega1)*cosd(i1)*cosd(theta1), ...
        -cosd(Omega1)*sind(i1);sind(i1)*sind(theta1), sind(i1)*cosd(theta1), cosd(i1)];

dv1_xyz = A_rth_xyz * dv1_rth;

%% 2

rplus_vec = rminus_vec;
rplus = norm(rplus_vec);
vplus_vec = vminus_vec + dv1_xyz;
vplus = norm(vplus_vec);

%% 3

gammaplus = asind(dot(rplus_vec,vplus_vec)/(rplus*vplus));

e2 = sqrt((((rplus*vplus^2/mu)-1)^2*cosd(gammaplus)^2)+sind(gammaplus)^2)

hhatplus = cross(rplus_vec,vplus_vec)/norm(cross(rplus_vec,vplus_vec));
rhatplus = rplus_vec/rplus;
thetahatplus = cross(hhatplus,rhatplus);

i2 = acosd(hhatplus(3));

Omega2a = asind(hhatplus(1)/sind(i2));
Omega2b = acosd(-hhatplus(2)/sind(i2));
Omega2 = Omega2a;

theta2a = asind(rhatplus(3)/sind(i2));
theta2c = 180 - asind(rhatplus(3)/sind(i2));
theta2b = acosd(thetahatplus(3)/sind(i2));
theta2d = -acosd(thetahatplus(3)/sind(i2));
theta2 = theta2c;

gammaplus = asind(dot(rplus_vec,vplus_vec)/(rplus*vplus))

thetastar2 = atan2d((rplus*vplus^2/mu*cosd(gammaplus)*sind(gammaplus)),(rplus*vplus^2/mu*(cosd(gammaplus))^2-1));

%% 4

de = e2-e1;
di = i2-i1;
dOmega = Omega2-Omega1;
dtheta = theta2-theta1;

%% 7

at1 = 6658;
at2 = 6798;
aiss = 6378+430;

TOF = pi*(sqrt(at1^3/mu)+sqrt(at2^3/mu));
niss = sqrt(mu/aiss^3);

phaseangle = 2*pi-niss*TOF;
phaseangle = phaseangle*180/pi;

