function f = f_equilibrium_2d(r,u,x,v,y,t,theta) 
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

f = 1./(exp( ((x-u).^2 + (y-v).^2)./t)./r + theta);
