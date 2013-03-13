function f = f_ES_equilibrium_1d(z,p,n,u,c,t,theta) 
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
% p: pressure
% n: density

f = 1./(exp( (c-u).^2 ./ t)./z + theta);
