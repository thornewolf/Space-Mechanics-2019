% Lagrange Coefficients
% Austin Kendall
% Last Updated: 19 FEB 2019
clear all
clc
%% Constants

mu = 42828; % km^2/s^2
a = 6482.3; % km
p = 1; % km

r1vec = [-2089.6 -2515.7 -6382]; % km
r1 = norm(r1vec); % km
r2 = 9416.57; % km

v1vec = [2.1744 0.76911 0.13452]; % km/s
v1 = norm(v1vec); % km/s

disp('Ran Constants')
%% Utilizing True Anomaly

thetastar1 = 360; % deg
thetastar2 = 67; % deg
dthetastar = thetastar2 - thetastar1; % deg

f = 1 - r1^2/p*(1 - cosd(dthetastar));

g = r2*r1/(sqrt(mu*p))*dthetastar; % km

fdot = dot(r1vec,v1vec)/(r1*p)*(1 - cosd(dthetastar)) - 1/r1*sqrt(mu/p)*sind(dthetastar); % 1/s

gdot = 1 - r1/p*(1 - cosd(dthetastar));

disp('Ran Utilizing True Anomaly')
%% Utilizing Eccentric Anomaly

E1 = -1.8076; % rad
E2 = 15.746; % rad
dE = E2 - E1; % rad
dt = 43200; % s

f = 1 - a/r1*(1 - cos(dE));

g = dt - sqrt(a^3/mu)*(dE - sin(dE)); % km

fdot = -sqrt(mu*a)/(r2*r1)*sin(dE); % 1/s

gdot = 1 - a/r2*(1 - cos(dE));

disp('Ran Utilizing Eccentric Anomaly')
%% Utilizing Hyperbolic Anomaly

H1 = 2*pi; % rad
H2 = 1.548764; % rad
dH = H2 - H1; % rad
dt = 2391; % s

f = 1 - abs(a)/r1*(cosh(dH) - 1);

g = dt - sqrt(abs(a)^3/mu)*(sinh(dH) - dH); % km

fdot = sqrt(mu*abs(a)/(r2*r1))*sinh(dH); % 1/s

gdot = 1 - abs(a)/r2*(cosh(dH) - 1);

disp('Ran Utilizing Hyperbolic Anomaly')
%% New Vector

r2vec = f*r1vec + g*v1vec; % km

v2vec = fdot*r1vec + gdot*v1vec; % km/s

fprintf('r2 = %d x %d y %d z km \n',r2vec)
fprintf('v2 = %d x %d y %d z km/s \n',v2vec)

disp('Ran New Vector')