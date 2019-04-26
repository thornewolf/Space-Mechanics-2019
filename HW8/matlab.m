clc; clear;
% Constants
planets = {'Sun', 'Moon', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptine', 'Pluto'};
rad_list = [695990.0, 1739.2, 2439.7, 6051.9, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 25269.0, 1162.0];
mu_list = [132712440000.0, 4902.8, 22032.0, 324860.0, 398600.0, 42828.0, 126713000.0, 37941000.0, 5794500.0, 6836500.0, 981.6];
sma_list = [0.0, 384400.0, 57910000.0, 108210000.0, 149600000.0, 227920000.0, 778570000.0, 1433530000.0, 2872460000.0, 4495060000.0, 5906380000.0];
R = containers.Map(planets,rad_list);
MU = containers.Map(planets,mu_list);
r = containers.Map(planets,sma_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
am = 28081;
AOPm = 229;
RAANm = 34.8;
em = 0.765;
incm = 20.6;
true_am = -45.6;

AOLm = AOPm + true_am;

vrm_eci = [-5978.4 -4668 -158.31]';
vvm_eci = [7.7465 -4.6198 -3.0877]';

dv = 4.5;
alpha = -30;
beta = 40;
pm = am*(1-em^2);
% FPA = acosd(norm(cross(vrm_eci, vvm_eci))/(norm(vrm_eci)*norm(vvm_eci)))
FPA = atan2d(norm(vrm_eci)*em*sind(true_am),pm);
phi = alpha + FPA;

% 1.

vdv_rth = dv*[cosd(beta)*sind(phi) cosd(beta)*cosd(phi) sind(beta)]';
vdv_eci = rot_rth_eci(RAANm, incm, AOLm) * vdv_rth;

% 2.
vrp_eci = vrm_eci;
vvp_eci = vvm_eci + vdv_eci;

% 3.
FPAp = asind(dot(vrp_eci,vvp_eci)/(norm(vrp_eci)*norm(vvp_eci)));
ep = sqrt((((norm(vrp_eci)*norm(vvp_eci)^2/MU('Earth'))-1)^2*cosd(FPAp)^2)+sind(FPAp)^2);

syms ap
eq = norm(vvp_eci)^2/2 - MU('Earth')/norm(vrp_eci) == MU('Earth')/(2*ap);
ap = double(solve(eq, ap));
pp = norm(cross(vrp_eci,vvp_eci))^2/MU('Earth');
h_hat = cross(vrp_eci, vvp_eci)/norm(cross(vrp_eci,vvp_eci));

incp = acosd(h_hat(3));

syms RAANp
eq(1) = h_hat(2) == -cosd(RAANp)*sind(incp);
eq(2) = h_hat(1) == sind(RAANp)*sind(incp);
RAANp1 = acosd(-h_hat(2)/sind(incp));
RAANp1 = [-RAANp1 RAANp1];
RAANp2 = asind(h_hat(1)/sind(incp));
RAANp2 = [180-RAANp2 RAANp2];

RAANp = min(RAANp2); %intersecting the two doesnt work due to small delta in solution. It's the negative value.

rp_eci_hat = vrp_eci/norm(vrp_eci);
theta_hat = cross(rp_eci_hat,h_hat);

AOLp1 = asind(rp_eci_hat(3)/sind(incp));
AOLp1 = [180-AOLp1 AOLp1];
AOLp2 = acosd(theta_hat(3)/sind(incp));
AOLp2 = [-AOLp2 AOLp2];
AOLp = min(AOLp2); %same as above. Small delta does not allow me to do intersection
AOLp = 191.1029;

true_ap = acosd((pp/norm(vrp_eci)-1)/ep); %Checked Gmat its neg
true_ap = 360 - true_ap;

% 4.
de = ep-em
dinc = incp-incm
dRAAN = RAANp-RAANm
dAOL = AOLp-AOLm

% 5.
vrp_eci; %-7192.67008574393           -3847.023676019902          355.568234258645            303.8898424306826  

% 6.
% I would not want to perform this maneuver because it results in a
% massively hyperbolic orbit. The goal of the maneuver is to do an orbital
% correction, not enter an escape trajectory. With an eccentricity of 2.02,
% the orbit is very far from bring elliptical. 

% 7.
at1 = 6658;
at2 = 6798;
a_sc = 6378+430;
TOF = pi*(sqrt(at1^3/MU('Earth'))+sqrt(at2^3/MU('Earth')));
n_sc = sqrt(MU('Earth')/a_sc^3);
phase = 2*pi-n_sc*TOF;
phase = phase*180/pi;

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