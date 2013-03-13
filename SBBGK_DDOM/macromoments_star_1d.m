function [n, nux, nE, ne] = macromoments_star_1d(J,k,w,f,v,ux)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(J.*w .* f);    % Density [n]
     nux = k*sum(J.*v .* w .* f);   % Macrospic moment [n*ux]
     nE  = k*sum(J.*( 0.5*v.^2 ).* w .* f);  % Energy Density [n*E]
     ne  = 0.5*k*sum(J.*( (v-ux).^2 ).* w .* f); % Internal energy [n*e]

