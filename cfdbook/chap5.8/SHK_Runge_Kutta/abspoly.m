% Amplification factor for SHK Runge-Kutta method
function res = abspoly(r)
global thetapar  center
% Coefficients for RK24 method Sommeijer et al.  Appl. Num. Meth. 16:201-225
c1 = 1/4; c2 = 1/3; c3 = 1/2; c4 = 1;
z = center + r*exp(i*thetapar);
res = - 1 + abs(1+c4*z*(1+c3*z*(1+c2*z*(1+c1*z))));

