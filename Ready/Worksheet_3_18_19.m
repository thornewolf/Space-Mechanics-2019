clc;
MU = 132712440000;

a_t = 791565000;
a_saturn = 1433530000;
mean_motion = sqrt(MU/a_saturn^3);


TOF = pi*sqrt(a_t^3/MU);
phase_angle = pi - mean_motion*TOF;
phase_angle_deg = rad2deg(phase_angle);

e_saturn = 0.0565;
a_t = 7.5107e8;
e_t = 0.80082;

syms E_arr
TOF = pi*sqrt(a_t^3/MU)
E_dep = -1.2449;
eq = sqrt(MU/a_saturn^3)*TOF == E_arr - e_saturn*sin(E_arr) - (E_dep - e_saturn*sin(E_dep));
E_arr = solve(eq, E_arr)
E_arr = 0;

true_a_dep = 2*atan(sqrt((1+e_saturn)/(1-e_saturn))*tan(E_dep/2))
true_a_arr = 2*atan(sqrt((1+e_saturn)/(1-e_saturn))*tan(E_arr/2))
dtrue_a = true_a_arr - true_a_dep
rad2deg(dtrue_a)
phase_angle = pi - dtrue_a
rad2deg(phase_angle)