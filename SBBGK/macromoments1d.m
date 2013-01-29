function [n, nux, E] = macromoments1d(k,w,f,v)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f);    % Density
     nux = k*sum(v .* w .* f);   % Density * velocity x
     E   = k*sum(1/2*( v.^2 ).* w .* f); % Total Energy Density