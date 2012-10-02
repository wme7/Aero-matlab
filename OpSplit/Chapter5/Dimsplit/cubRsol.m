function [u,s] = cubRsol(ul,ur,delta)
%
% A wrapper for a scalar Riemannsolver
% to fit with genFrontTrack(...). Uses the flux function 
% f=scalflux(u);
%
flux='cub';

[u,s]=ScalarRiemannsol(ul,ur,delta,flux);


end

