% AE 313 Hw 6 Ground Track
clear all; clc;
addpath ../Functions
Re = 6378; % [km] radius of Earth
mu = 398600; % [km^3/s^2] grav param of Earth
a = 6*Re; % [km] semi-major axis
e = 0.7; % eccentricity
w = 25*pi/180;% [rad] argument of periapsis
Om = 220*pi/180; % [rad] longitude of ascending node
i = 30*pi/180; % [rad] inclination
th = 10*pi/180;% [rad] argument of latitude
JD = 2457997.5; % [days] Julian Date
ts1 = 290*pi/180; % [rad] new true anomaly
Tu = JD - 2451545;
ts = th - w; % [rad] previous true anomaly
M = E - e*sin(E);
t_0p = sqrt(a^3/mu)*M; % [sec] previous time since periapsis
E1 = 2*atan(sqrt((1-e)/(1+e))*tan(ts1/2)); % [rad] new Eccentric anomlay
M1 = E1 - e*sin(E1);
t_1p = sqrt(a^3/mu)*M1; % [sec] time since periapsis
dts = ts1 - ts;
dt = t_1p - t_0p + Period; % [sec] time difference between the two
if ts1 > 0 && ts > ts1
 dts = 2*pi - ts + ts1; % [rad] change in true anomaly
end
tS = [ts:.1.*pi/180:ts1]; % [rad] true anomaly array
% Eccentric anomaly must be a smooth line
EV = mod(2.*atan(sqrt((1-e)/(1+e)).*tan(tS/2)),2*pi); % [rad] ecc. anom
t_V0 = sqrt(a^3/mu).*(EV - e.*sin(EV) - E + e*sin(E)); % [sec] time array
TuV = Tu + t_V0./3600./24; % [day] J200 date
thERAV = 2*pi*(0.7790572732640 + 1.00273781191135448.*TuV);
% Latitude and Longitude
LatV = asin(sin(i).*sin(tS + w));
alpha = atan2(cos(i).*sin(tS + w), cos(tS + w));
LongV = alpha + Om - thERAV;
% Plot results
figure; hold on; set(gca,'fontsize',14)
 geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
 geoshow(LatV.*180./pi,LongV.*180./pi)
 xlabel('Longitude \lambda, deg'); ylabel('Latitude \phi, deg')