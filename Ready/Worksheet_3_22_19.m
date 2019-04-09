clc;
MU = 398600;
r_earth = 6378;

am = 4000 + r_earth;
em = 0.3;
tam = 50;

ap = 20000;
ep = 0.4;
tap = 200;

syms rm

pm = am*(1-em^2);
rm = pm/(1+em*cosd(tam))

syms FPAm
eq = tand(FPAm) == rm*em/pm*sind(tam);
FPAm = double(solve(eq, FPAm))

hm = sqrt(MU*pm)

syms vm
eq = hm == rm*vm*cosd(FPAm);
vm = double(solve(eq, vm))









