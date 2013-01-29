function [n, nux, E] = macromoments_star_1d(J,k,w,f,v)
%% Compute Macroscopic Moments
% Using Quadrature rules to integrate for:
     n   = k*sum(J.*w .* f);    % Density [n]
     nux = k*sum(J.*v .* w .* f);   % Macrospic moment [n*ux]
     E   = k*sum(J.*( 1/2*v.^2 ).* w .* f); % Energy Density [E]

