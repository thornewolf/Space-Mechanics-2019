syms r_1 vr_1 vr_2 r_2 v_1 v_2 vv_1 vv_2 P true_2 true_1 mu a p E
vr_1 = sym('vr',[1 2])
vv_1 = sym('vv',[1 2])
f = 1-r_2/p*(1-cosd(true_2-true_1)); 
g = r_2*r_1/sqrt(mu*p)*sind(true_2-true_1);
f_dot_general = (vr_1(1)*vv_1(1)+vr_1(2)+vv_1(2))%/(p*r_1)*(1-cosd(true_2-true_1)) - 1/r_1*(mu/p)^0.5*sind(true_2-true_1)
% f_dot_elliptic = -(mu*a)^0.5/(r_1*r_2)*sind(E_2-E_1);
g_dot = 1 - r_1/p*(1-cosd(true_2-true_1));

    f_dot_sol = subs(f_dot_general,...
{mu p r_1 vr_1(1) vr_1(2) vv_1(1) vv_1(2) true_2 true_1},...
{398600 51024 25119.5 12559.8, 21754.1 -3.0257, 3.8431 190 60})

double(f_dot_sol)