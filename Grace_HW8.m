%% AE 313 HW 8
% Grace Day
clc; close all
a = 28081; %km
w = 229; %deg
Omega = 34.8; %deg
e = 0.765;
i = 20.6; %deg
TS = -45.6; %deg
rminus = [-5978.4 -4668 -158.31]; %xyz km
vminus = [7.7465 -4.6198 -3.0877]; %xyz km/s
deltav = 4.5; %km/s
alpha = -30; %deg
beta = 40; %deg
mu = 398600;
% Problem 1 
    % Find deltav in rOh and xyz
p = a*(1-(e^2));
gamma1 = atan2d((((norm(rminus)*e)/p)*sind(TS)),1);
phi = alpha + gamma1;
dVROH = deltav*[(cosd(beta)*sind(phi)) (cosd(beta)*cosd(phi)) ...
    (sind(beta))];
theta = TS + w;
dVXYZ =(rot_xyz_rth(Omega,i,theta)).*dVROH;
dVXYZ = [dVXYZ(1)+dVXYZ(4)+dVXYZ(7) dVXYZ(2)+dVXYZ(5)+dVXYZ(8)...
    dVXYZ(3)+dVXYZ(6)+dVXYZ(9)]
% Problem 2
    % find rplus and vplus in xyz
rplus = rminus;
vplus = vminus + dVXYZ;
% Problem 3 
    % Find eplus, iplus, OmegaPlus, thetaplus, and TAplus
hplus = (cross(rplus,vplus))/(norm(cross(rplus,vplus)));
iplus = acosd(hplus(3));
OmegaplusS = asind((hplus(1)/sind(iplus)));
OmegaplusC = -acosd((hplus(2))/-sind(iplus)); %right one
rhatplus = rplus/norm(rplus);
thetahatplus = cross(hplus,rhatplus);
thetaplusS = asind(rhatplus(3)/sind(iplus)); %right one
thetaplusC = acosd(thetahatplus(3)/sind(iplus));
gammaplus = asind((dot(rplus,vplus))/(norm(rplus).*norm(vplus)));
    %less than 0 so descending
TSplus = atan2d(((norm(rplus)*(norm(vplus)^2))/mu)*cosd(gammaplus)*...
    sind(gammaplus),(((norm(rplus)*(norm(vplus)^2))/mu)*...
    cosd(gammaplus)^2)-1);
eplus = sqrt((((norm(rplus)*(norm(vplus)^2))/mu)-1)^2*...
    (cosd(gammaplus)^2)+(sind(gammaplus)^2));
% Problem 4
    % Find deltae, deltai, deltaOmega, and deltatheta
deltae = eplus-e;
deltai = iplus-i; %deg
deltaOmega = OmegaplusC-Omega+360;
deltatheta = thetaplusS-theta;
% Problem 7
aT1 = 6658; %km
aT2 = 6798; %km
TA = 2*pi;
mu = 398600;
n = sqrt(mu/(aT2^3)); %rad/sec
TOF = pi*sqrt((aT1^3)/mu) + pi*sqrt((aT2)/mu); %sec
PA = TA - (n*TOF); %deg

function A = rot_xyz_rth(o,i,t)
A = inv(rot_rth_eci(o,i,t));
end

function A = rot_rth_eci(o,i,t)

% o : Omega, Longitude of the Ascending Node (RAAN)
% i : i, inclination
% t : theta, Argument of Latitude

A = [cosd(o)*cosd(t)-sind(o)*cosd(i)*sind(t), -cosd(o)*sind(t)-sind(o)*cosd(i)*cosd(t), sind(o)*sind(i); ...
        sind(o)*cosd(t)+cosd(o)*cosd(i)*sind(t), -sind(o)*sind(t)+cosd(o)*cosd(i)*cosd(t), -cosd(o)*sind(i); ...
        sind(i)*sind(t), sind(i)*cosd(t), cosd(i)];

end