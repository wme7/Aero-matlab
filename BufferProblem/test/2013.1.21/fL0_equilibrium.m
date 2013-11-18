function Ml0 = fL0_equilibrium(rl0,ul0,c,tl0,theta) 
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

 Ml0 = 1./(exp( (c-ul0).^2 ./ tl0)./rl0 + theta);
 
