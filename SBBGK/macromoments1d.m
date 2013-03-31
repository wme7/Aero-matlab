function [n, nux, nE ,ne , Wxx] = macromoments1d(k,w,f,v,ux)
%% Compute Macroscopic Moments

% For any finite difference method (FDM), we use normal sum().
% Using Quadrature rules to integrate for:
     n   = k*sum(w .* f );    % Density [n]
     nux = k*sum(v .* w .* f );   % momentum in x [n*ux]
     nE  = 0.5*k*sum(( v.^2 ).* w .* f ); % Total Energy [n*E]
     
	 %ux = nux/n; <-- Carefull! must use 'ux' from the last time step!
	 ne  = 0.5*k*sum(( (v-ux).^2 ).* w .* f); % Internal Energy [n*e]
     Wxx = k*sum(w .* f .*(v-ux).^2 ); % Pressure tensor component xx
     

