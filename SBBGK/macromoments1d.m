function [n, j_x, E] = macromoments1d(k,w,f,v)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f);    % Number density
     j_x = k*sum(v .* w .* f);   % Macrospic moment in x
     E   = k*sum(1/2*( v.^2 ).* w .* f); % Total Energy Density