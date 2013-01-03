% Amplification factor for Runge-Kutta method of Le and Moin JCP92:369-379
function res = abspoly(r)
global thetapar  center
% Coefficients for Runge-Kutta method of Le and Moin: 
sig1 = 8/15; sig2 = 5/12; sig3 = 3/4;
zeta2 = -17/60; zeta3 = -5/12; 
alp1 = 4/15; alp2= 1/15; alp3 = 1/6;
bet1 = 4/15; bet2= 1/15; bet3 = 1/6;
z= center + r*exp(i*thetapar);
v = real(z); w = imag(z);
R1 = (1+alp1*v+i*sig1*w)/(1-bet1*v);
R2 = (1+alp2*v+i*sig2*w)/(1-bet2*v);
R3 = (1+alp3*v+i*sig3*w)/(1-bet3*v);
pol = R3*(R2*R1 + i*zeta2*w/(1-bet2*v)) + i*zeta3*w*R1/(1-bet3*v);
res = - 1 + abs(pol);
