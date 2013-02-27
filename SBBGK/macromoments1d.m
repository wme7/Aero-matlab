function [n, nux, E] = macromoments1d(k,w,f,v)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f ,3);    % Density
     nux = k*sum(v .* w .* f ,3);   % Density * velocity x
     E   = k*sum(0.5*( v.^2 ).* w .* f ,3); % Total Energy Density