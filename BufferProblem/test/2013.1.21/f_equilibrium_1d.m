function f = f_equilibrium_1d(r,u,c,t,theta) 
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

 f = 1./(exp( (c-u).^2 ./ t)./r + theta);
 
% r.*ppp
% f = r.*(1/2*pi*RRR*t)./(exp(- ((x-u).^2) ./ (2*RRR*t)))  