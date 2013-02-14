function f = f_equilibrium_2d(r,ux,uy,t,vx,vy,theta) 
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
f = 1./(((1./r).* exp( ((vx-ux).^2 + (vy-uy).^2 )./t) + theta));


