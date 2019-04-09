%% AE 313 Group Exam 2
% Team 3 - Jiyoung Hwang, Grace Day, Thorne Wolfenbarger, Aaron Scott

clc; clear all; close all;

% Knowns
MU = 398600;
a = 42170; %km
e = 0.01053;
RAAN = 14.89; %deg
AOP = 318.4; %deg
inc = 14.37; %deg
true_a1 = 184.2; %deg
energy = -MU/(2*a);

% 1. Find the current position and velocity in ECI and perifocal coordinates
p = a*(1-e^2);
h = sqrt(MU*p);
r1 = p/(1+e*cosd(true_a1));

vr1_rth = [r1 0 0]';
vv1_rth = [h*e/p*sind(true_a1) h/r1 0]';

AOL1 = AOP + true_a1;
rth_eci = rot_rth_eci(RAAN,inc,AOL1);

vr1_eci = rth_eci*vr1_rth;
vv1_eci = rth_eci*vv1_rth;

vr1_peri = r1*[cosd(true_a1) sind(true_a1) 0]';
vv1_peri = [-h/p*sind(true_a1) h/p*(e+cosd(true_a1)) 0]';
% 2. You need to send the message by the time E=15deg. How long do you have
% until t2?
period = 2*pi*sqrt(a^3/MU);
E2 = 15; %deg
E2_rad = E2*pi/180;

% syms E1
% eq = r1 == a*(1-e*cosd(E1));
% E1 = double(solve(eq,E1))
E1 = 2*atand(tand(true_a1/2)/sqrt((1+e)/(1-e)));
% E1 = acosd((a*e+r1*cosd(true_a1))/a);
% E1 = acosd(-(r1/a-1)/e)*sign(true_a1);
E1_rad = E1*pi/180;

dt = sqrt(a^3/MU)*(E2_rad - E1_rad - (e*sin(E2_rad) - e*sin(E1_rad)));

% 3. What is the new position and velocity in ECI and perifocal coordinates?
r2 = a*(1-e*cos(E2_rad));
true_a2 = 2*atan2d(sqrt((1+e)/(1-e))*tan(E2_rad/2),1);

% f = 1 - a/r1*cos(E2_rad - E1_rad);
% g = dt - sqrt(a^3/MU)*((E2_rad - E1_rad) - sin(E2_rad - E1_rad));
% fdot = sqrt(MU*a)/(r2*r1)*sin(E2_rad - E1_rad);
% gdot = 1-a/r2*(1-cos(E2_rad - E1_rad));
% vr2_eci = f*vr1_eci + g*vv1_eci;
% vv2_eci = fdot*vr1_eci + gdot*vv1_eci;
vr2_rth = [r2 0 0]';
vv2_rth = [h*e/p*sind(true_a2) h/r2 0]';

AOL2 = AOP + true_a2;
rth_eci = rot_rth_eci(RAAN,inc,AOL2);

vr2_eci = rth_eci*vr2_rth;
vv2_eci = rth_eci*vv2_rth;

vr2_peri = r2*[cosd(true_a2) sind(true_a2) 0]';
vv2_peri = [-h/p*sind(true_a2) h/p*(e+cosd(true_a2)) 0]';


% 4. Find the latitiude and longitide at the two positions
dt_days = dt/60/60/24;
rot_earth = 7.2921151467*10^-5; %rad
tu1 = juliandate(datetime('2019-03-06 03:00:00')) - 2451545;
tu2 = juliandate(datetime('2019-03-06 03:00:00')) - 2451545;

theta_era1 = 2*pi*(0.7790572732640 + 1.00273781191135448*tu1);
theta_era1 = mod(theta_era1,2*pi)*180/pi;
theta_era2 = 2*pi*(0.7790572732640 + 1.00273781191135448*tu2);
theta_era2 = mod(theta_era2,2*pi)*180/pi;

eci_ecef1 = rot_eci_ecef(theta_era1);
eci_ecef2 = rot_eci_ecef(theta_era2);
vr1_ecef = eci_ecef1*vr1_eci;
vr2_ecef = eci_ecef2*vr2_eci;


alpha1 = atan2d(cosd(inc)*sind(AOL1),cosd(AOL1));
alpha2 = atan2d(cosd(inc)*sind(AOL2),cosd(AOL2));
theta_gr1 = theta_era1+180/pi*rot_earth*0;
theta_gr2 = theta_era2+180/pi*rot_earth*dt;

lat1 = asind(vr1_ecef(3)/r1);
long1 = alpha1 + RAAN - theta_gr1;
lat2 = asind(vr2_ecef(3)/r2);
long2 = alpha2 + RAAN - theta_gr2;
long2 = mod(long2,-360);

% 5. Create a ground track between the two positions using MATLAB. Mark the
% starting location.

true_a_array = linspace(true_a1, true_a2+360, 1000);
r_array = p./(1+e*cosd(true_a_array));
% E_array = acosd(-(r_array/a-1)/e).*sign(true_a_array);
E_array = 2*atand(tand(true_a_array/2)/sqrt((1+e)/(1-e)));
E_array_rad = E_array*pi/180;

dt_array = sqrt(a^3/MU)*(E_array_rad - E1_rad - (e*sin(E_array_rad) - e*sin(E1_rad)));

theta_era = theta_era1;
theta_gr_array = theta_era+180/pi*rot_earth*dt_array;


theta_array = AOP + true_a_array;
lat_array = asind(sind(inc)*sind(theta_array));

alpha_array = atan2d(cosd(inc)*sind(theta_array),cosd(theta_array));
long_array = alpha_array + RAAN - theta_gr_array;

figure; hold on;
geoshow("landareas.shp","FaceColor",[0.5 1.0 0.5]);
geoshow(lat_array, long_array,'Color','red');
geoshow(lat_array(1),long_array(1),'DisplayType','Point','Marker','x','Markersize',20);
xlabel('long (\circ)')
ylabel('lat (\circ)')

% 6. Assume you have the ground station in Prescott. Plot the elevation angle
% versus time passed. If the ground station needs a minimum elevation angle
% of 10 degrees, how much time is the TDRS visible? Can you send your
% urgent message to the TDRS?
min_elevation = 10; %deg
radius_e = 6378; %km
lat_gs = 34.54; %deg
long_gs = 112.4685; %deg
h_gs = 1.64; %km
vr_e_gs = [0 0 radius_e+h_gs]'; %km

elevation_array = zeros(1,length(r_array));
for i = 1:length(r_array)
    vr_rth = [r_array(i) 0 0]';
    rth_eci = rot_rth_eci(RAAN, inc, theta_array(i));
    eci_ecef = rot_eci_ecef(theta_era);
    vr_ecef = eci_ecef*rth_eci*vr_rth;

    ecef_sez = rot_ecef_sez(long_gs, lat_gs);
    vr_sez = ecef_sez*vr_ecef;

    vr_e_sc = vr_sez;

    vr_gs_sc = vr_e_sc - vr_e_gs;

    elevation_array(i) = asind(vr_gs_sc(3)/norm(vr_gs_sc));
end

figure(2)
hold on;
plot(dt_array, elevation_array)
line([0 dt_array(end)], [10 10],'Color',[0.78 0.43 0.02])
xlabel('Time (s)')
ylabel('Elevation (\circ)')
axis([0 dt_array(end) elevation_array(1) 40])
q = 1;
while elevation_array(q)<10
    q = q + 1;
end
time_start = dt_array(q);
time_end = dt_array(end);
time_above_minimum = time_end - time_start;
time_above_minimum_hours = time_above_minimum/60/60;
% 7. Why would the TDRS be in such an orbit?



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
