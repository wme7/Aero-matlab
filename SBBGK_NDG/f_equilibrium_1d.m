function f = f_equilibrium_1d(z,u,c,t,theta) 
% Compute equilibrum in 1d cases for:
%
% MB:  Maxwell-Boltzmann, theta =  0
% FD:  Fermi-Diract,      theta = +1
% BE:  Bose-Einstein,     theta = -1
%
% Inputs:
% u: macroscopic velocity
% x: microscopic velocity
% t: temperature
% z: fugacity

f = 1./(exp( (c-u).^2 ./ t)./z + theta);
