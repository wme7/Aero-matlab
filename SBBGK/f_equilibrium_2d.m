function f = f_equilibrium_2d(r,u1,u2,v1,v2,t,theta) 
% Compute equilibrum in 1d cases for:
%
% MB:  Maxwell-Boltzmann, theta =  0
% FD:  Fermi-Diract,      theta = +1
% BE:  Bose-Einstein,     theta = -1
%
% inputs:
% u: macroscopic velocity
% x: microscopic velocity
% v: macroscopic velocity
% y: microscopic velocity
% t: temperature
% r: fugacity
%
f = 1./((1./r).*(exp( ((v1-u1).^2 + (v2-u2).^2 ) ./t) + theta));

