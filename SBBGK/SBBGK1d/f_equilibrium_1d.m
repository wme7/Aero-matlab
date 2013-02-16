function f = f_equilibrium_1d(r,u,x,t,theta) 
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

f = 1./(exp( (x-u).^2 ./ t)./r + theta);
