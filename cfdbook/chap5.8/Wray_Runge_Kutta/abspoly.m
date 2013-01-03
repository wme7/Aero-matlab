% Amplification factor for Wray's Runge-Kutta method
function res = abspoly(r)
global thetapar  center
% Coefficients for Wray's Runge-Kutta method: Spalart et al. JCP 96:297 1991
%gamma1 = 8/15; gamma2 = 5/12; gamma3 = 0.75;
%zeta1 = -17/60; zeta2 = -5/12; 
z = center + r*exp(i*thetapar);
%pol = (1+gamma3*z)*(1+gamma2*z)*(1+gamma1*z);
%pol = pol + zeta1*z*(1+gamma3*z) + zeta2*z*(1+gamma1*z);
pol = 1 + z*(1 + z*(0.5 + z/6));
res = - 1 + abs(pol);
