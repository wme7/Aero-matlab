function [un,up] = WENO3_1D_flux(u)
% This function assumes that u = [u(1) u(2) u(3) u(4) u(5)]
% And we return: up: u_i+1/2^(+) & un: u_i+1/2^(-)

% Constants
dn0 = 1/10; dn1 = 6/10; dn2 = 3/10; dp0 = 3/10; dp1 = 6/10; dp2 = 1/10;
epsilon = 1E-7;

%% Compute WENO3 1d Fluxes:
% ENO polynomials
pn0  = (-u(1) +  2*u(2) + 23*u(3)) / 24;
pn1  = ( u(1) -  4*u(2) +  3*u(3)) / 2 ;
pn2  = ( u(1) -  2*u(2) +    u(3)) / 2 ;

pm0  = (-u(2) + 26*u(3) -    u(4)) / 24;
pm1  = (-u(2)           +    u(4)) / 2 ;
pm2  = ( u(2) -  2*u(3) +    u(4)) / 2 ;

pp0  = (-u(3) +  2*u(4) + 23*u(5)) / 24;
pp1  = (-u(3) +  4*u(4) -  3*u(5)) / 2 ;
pp2  = ( u(3) -  2*u(4) +    u(5)) / 2 ;

% Nonlinear Smoothness Indicators
ISn = 13/3*pn2^2 + pn1^2;
ISm = 13/3*pm2^2 + pm1^2;
ISp = 13/3*pp2^2 + pp1^2;

%% Positive Flux
% Non-normalized stencil weights
alphapn = dp0 /(epsilon + ISn)^2;
alphapm = dp1 /(epsilon + ISm)^2;
alphapp = dp2 /(epsilon + ISp)^2;
alphasump = alphapn + alphapm + alphapp;

% Weigths
wpn = alphapn / alphasump;
wpm = alphapm / alphasump;
wpp = alphapp / alphasump;

% Negative Flux
% Non-normilized Stencil weights
alphann = dn0 /(epsilon + ISn)^2;
alphanm = dn1 /(epsilon + ISm)^2;
alphanp = dn2 /(epsilon + ISp)^2;
alphasumn = alphann + alphanm + alphanp;

% Weigths
wnn = alphann / alphasumn;
wnm = alphanm / alphasumn;
wnp = alphanp / alphasumn;

%% WENO Fluxes
% Positive Flux: u_i+1/2 (+)
dx  =-0.5;
up  = wpn * (pn0 + pn1*dx + pn2*dx^2) ...
	+ wpm * (pm0 + pm1*dx + pm2*dx^2) ...
	+ wpp * (pp0 + pp1*dx + pp2*dx^2);

% Negative Flux: u_i+1/2 (-)
dx  = 0.5;
un  = wnn * (pn0+ pn1*dx + pn2*dx^2) ...
    + wnm * (pm0+ pm1*dx + pm2*dx^2) ...
	+ wnp * (pp0+ pp1*dx + pp2*dx^2);

end