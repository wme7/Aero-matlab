function f = f_ES_equilibrium_1d(r,p,n,u,x,t,theta) 
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

f = 1./((1./r).*exp((x-u).^2./t) + theta);
b = p/n;
