function f = f_equilibrium_star_1d(r,c,theta) 
% Compute equilibrum in 1d cases for:
%
% MB:  Maxwell-Boltzmann, theta =  0
% FD:  Fermi-Diract,      theta = +1
% BE:  Bose-Einstein,     theta = -1
%
% inputs:
% u: macroscopic velocity
% x: microscopic velocity
% t: temperature
% r: fugacity

f = 1./(r.*exp(c.^2) + theta);
