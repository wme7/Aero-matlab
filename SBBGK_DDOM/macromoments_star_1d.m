function [n, nux, nE, ne, Wxx] = macromoments_star_1d(J,k,w,f,v,ux)
%% Compute Macroscopic Moments

% For any finite difference method (FDM), we use normal sum().
% Using Quadrature rules to integrate for:
     n   = k*sum(J.*w .* f);    % Density [n]
     nux = k*sum(J.*v .* w .* f);   % Macrospic moment [n*ux]
     nE  = k*sum(J.*( 0.5*v.^2 ).* w .* f);  % Energy Density [n*E]
     %ux = nux/n; <-- Carefull! must use 'ux' from the last time step!
     ne  = 0.5*k*sum(J.*( (v-ux).^2 ).* w .* f); % Internal energy [n*e]
     Wxx = k*sum(J.*w .* f .*(v-ux).^2 ); % Pressure tensor component xx

