% Specification of Riemann problems

global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  lambda

%  Shocktube problem of G.A. Sod, JCP 27:1, 1978 
pleft = 1.0;  pright = 0.1; rholeft = 1.0;  rhoright = 0.125;
uleft = 0;  uright = 0; tend = 0.15; lambda = 0.5;	% lambda = dt/dx

% Lax test case: M. Arora and P.L. Roe: JCP 132:3-11,  1997
% pleft = 3.528;  pright = 0.571; rholeft = 0.445;  rhoright = 0.5;
% uleft = 0.698;  uright = 0; tend = 0.15; lambda = 0.2;	% lambda = dt/dx

% Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11,  1997
% pleft = 10.333;  pright = 1; rholeft = 3.857;  rhoright = 1;
% uleft = 0.92;  uright = 3.55; tend = 0.09; lambda = 0.2; % lambda = dt/dx

% Shocktube problem with supersonic zone
%pleft = 1;  pright = 0.02; rholeft = 1;  rhoright = 0.02;
%uleft = 0;  uright = 0; tend = 0.162; lambda = 0.3;	% lambda = dt/dx 

% Contact discontinuity
%pleft = 0.5;  pright = 0.5; rholeft = 1.0;  rhoright = 0.6;
%uleft = 0;  uright = 0; tend = 1; lambda = 0.4; 	% lambda = dt/dx
