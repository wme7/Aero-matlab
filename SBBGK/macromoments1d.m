function [n, nux, E, Wxx] = macromoments1d(k,w,f,v,ux)
%% Compute Macroscopic Moments

% For any finite difference method (FDM), we use normal sum().
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f );    % Density
     nux = k*sum(v .* w .* f );   % Density * velocity x
     E   = k*sum(0.5*( v.^2 ).* w .* f ); % Total Energy Density
     %ux = nux/n; <-- Carefull! must use 'ux' from the last time step!
     Wxx = k*sum(w .* f .*(v-ux).^2 ); % Pressure tensor component xx

