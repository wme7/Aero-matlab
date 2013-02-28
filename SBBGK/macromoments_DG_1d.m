function [n, nux, E] = macromoments_DG_1d(k,w,f,v)
%% Compute Macroscopic Moments

% For DG Method (sum in the 3rd direction of our arrays)
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f ,3);    % Density
     nux = k*sum(v .* w .* f ,3);   % Density * velocity x
     E   = k*sum(0.5*( v.^2 ).* w .* f ,3); % Total Energy Density
     %Wxx = k*sum(w .* f .*(v-u).^2 ,3); % Pressure tensor component xx