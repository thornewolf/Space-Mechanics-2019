%% AE 313 Group Exam 1
% Team 3 - Jiyoung Hwang, Grace Day, Thorne Wolfenbarger, Aaron Scott

% Position Vectors
rSunMars = [198720000, 72978000, -334700]; % km
rMarsPhobos = [174.29, -9319.9, -702.19]; % km
rMarsInSight = [-236120, 638170, -17906]; % km
rSunPhobos = rSunMars + rMarsPhobos; % km
rSunInSight = rSunMars + rMarsInSight; % km
rPhobosInSight = rMarsInSight - rMarsPhobos; % km
rInSightSun = -rSunInSight; % km
rMarsSun = -rSunMars; % km
rInSightPhobos = -rPhobosInSight; % km

% Gravitational Parameters
muInSight = 4.6318e-17; % km^3/s^2
muPhobos = 0.000709; % km^3/s^2
muSun = 132712440000; % km^3/s^2
muMars = 42828; % km^3/s^2
G = 6.67*10^-11; % N*m^2/kg^2

% Masses
M_InSight = (muInSight) / G; % kg
M_Phobos = (muPhobos) / G; % kg
M_Sun = (muSun) / G; % kg
M_Mars = (muMars) / G; % kg

%% Q1 - 
% Calculate the center of mass with respect to the Sun

rCM = (rSunMars * M_Mars + rSunPhobos * M_Phobos + rSunInSight * M_InSight) ...
      / (M_InSight + M_Phobos + M_Mars + M_Sun);

%% Q2 -
% Find the acceleration of InSight. What is the acceleration relative to?

Accel_InSight = -G * (((M_Sun * rSunInSight) / norm(rSunInSight)^3 ) + ...
                     ((M_Phobos * rPhobosInSight) / (norm(rPhobosInSight)^3 )) + ...
                     ((M_Mars * rMarsInSight) / (norm(rMarsInSight)^3))); 

% The Acceleration is relative to the center of mass (inertial point) 


%% Q3 - 
% Write the differential equation that governs the relative motion of
% InSight relative to Mars.

% Answers in Document

%% Q4 - 
% Determine the dominant, direct, and indirect acceleration (vector)
% on InSight
fprintf(['Question 4'])
Accel_InSight_due_to_Mars = (-G * (M_InSight + M_Mars) ...
                            / (norm(rMarsInSight)^3) * (rMarsInSight))

Dir_Accel_InSight_due_to_Sun = (G * (M_Sun * ((rInSightSun / norm(rInSightSun) ^3))))
InDir_Accel_InSight_due_to_Sun = (G * (M_Sun * (-(rMarsSun / norm(rMarsSun)^3))))

Dir_Accel_InSight_due_to_Phobos = (G * (M_Phobos * ((rInSightPhobos / norm( rInSightPhobos)^3))))
InDir_Accel_InSight_due_to_Phobos = (G * (M_Phobos * (-(rMarsPhobos/norm(rMarsPhobos)^3))))
% Answers in Document

%% Q5 - 
% For the same relative acceleration (InSight relative to Mars), evaluate
% the acceleration due to Mars, acceleration due to the Sun, and
% acceleration due to Phobos

Accel_InSight_due_to_Mars = (-G * (M_InSight + M_Mars) ...
                            / (norm(rMarsInSight)^3) * (rMarsInSight));
                        
Accel_InSight_due_to_Sun = (G * (M_Sun * ((rInSightSun / norm(rInSightSun)^3) ...
                           - (rMarsSun / norm(rMarsSun)^3))));
                        
Accel_InSight_due_to_Phobos = (G * (M_Phobos * ((rInSightPhobos / norm(rInSightPhobos)^3) ...
                              - (rMarsPhobos/norm(rMarsPhobos)^3))));

%% Q6 - 
% At this instant, is it reasonable to assume a two-body relative motion
% for motion of the spacecraft with respect to Mars? Why or why not?

% It is not reasonable to assume a 2-body relative model for motion
% of the spacecraft with respect to Mars because the gravitational
% influence of the Sun is only 1 order of magnitude less than that of Mars.
% This is a 10% error, which can lead to an unpredictable trajectory and
% landing zone."

fprintf(['\nQuestion 1: rCM = <%.3e, %.3e, %.3e> km\n' ...
         'Question 2: a_InSight = <%.3e, %.3e, %.3e> km/s^2\n' ...
         'Question 5: a_InSight_Due_To_Mars = <%.3e, %.3e, %.3e> km/s^2\n' ...
         '            a_InSight_Due_To_Sun = <%.3e, %.3e, %.3e> km/s^2\n' ...
         '            a_InSight_Due_To_Phobos = <%.3e, %.3e, %.3e> km/s^2\n'], ...
         rCM, Accel_InSight, Accel_InSight_due_to_Mars, ...
         Accel_InSight_due_to_Sun, Accel_InSight_due_to_Phobos);