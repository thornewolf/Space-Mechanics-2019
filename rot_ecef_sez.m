function A = rot_ecef_sez(l,p)

% l : lamda_gs, longitude of the ground station
% p : phi_gs, latitude of the ground station

A = [sind(p)*cosd(l), sind(p)*sind(l), -cosd(p);
    -sind(l), cosd(l), 0;
    cosd(p)*cosd(l), cosd(p)*sind(l), sind(p)];
    
end