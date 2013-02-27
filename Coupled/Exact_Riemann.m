function [xx,rhoexact,uexact,pexact,machexact,entroexact,energexact] = ...
    Exact_Riemann(rholeft,uleft,pleft,rhoright,uright,pright,tend)
%---------------------------------------------
% Riemann Solver for solving shoc-tube problem
% by Manuel Diaz, 03.09.2011
%---------------------------------------------
% Institute of Applied Mechanics
% Aerodynamic Design and Analsis Lab
% Laboratory 007
%---------------------------------------------
% NOTE:
% A Cavitation Check is the is incorporated in
% the code. It further prevents plotting for
% possible but physically unlikely case of ex-
% Expansion shocks.
%---------------------------------------------
% INPUT VARIABLES:
% Problem definition: Conditions at time t=0
% u1, p1, rho1
% u4, p4, rho4
global PRL  CRL MACHLEFT  gamma

gamma = 2.5; gammab = 1/(gamma - 1); gam1 = gamma-1;

% Assumed structure of exact solution
%
%    \         /      |con |       |s|
%     \   f   /       |tact|       |h|
% left \  a  /  state |disc| state |o| right
% state \ n /    2    |cont|   3   |c| state
%   1    \ /          |tinu|       |k|   4
%         |           |ity |       | |

PRL = pright/pleft;
cright = sqrt(gamma*pright/rhoright); cleft = sqrt(gamma*pleft/rholeft);
CRL = cright/cleft;
MACHLEFT = (uleft - uright)/cleft;

p34 = fzero('f',3);		% p34 = p3/p4
p3 = p34*pright; 	alpha = (gamma+1)/(gamma-1);
rho3 = rhoright*(1+alpha*p34)/(alpha+p34); 
rho2 = rholeft*(p34*pright/pleft)^(1/gamma);
u2 = uleft-uright+(2/(gamma-1))*cleft*...
	(1-(p34*pright/pleft)^((gamma-1)/(2*gamma)));
c2 = sqrt(gamma*p3/rho2);
spos = 0.5 + ...		% Shock position
	tend*cright*sqrt((gamma-1)/(2*gamma) + (gamma+1)/(2*gamma)*p34)+...
	tend*uright;

conpos = 0.5 + u2*tend + tend*uright;	% Position of contact discontinuity 
pos1 = 0.5 + (uleft - cleft)*tend;	% Start of expansion fan
pos2 = 0.5 + (u2+uright-c2)*tend;	% End of expansion fan
xx = 0:0.002:1;
pexact = zeros(size(xx)); uexact= zeros(size(xx)); rhoexact = zeros(size(xx));
machexact = zeros(size(xx));  cexact = zeros(size(xx));
for i = 1:length(xx)
  if xx(i) <= pos1
    pexact(i) = pleft;    rhoexact(i) = rholeft;
    uexact(i) = uleft;    cexact(i)   = sqrt(gamma*pexact(i)/rhoexact(i));
    machexact(i) = uexact(i)/cexact(i);
  elseif xx(i) <= pos2
    pexact(i) = pleft*(1+(pos1-xx(i))/(cleft*alpha*tend))^(2*gamma/(gamma-1));
    rhoexact(i) = rholeft*(1+(pos1-xx(i))/(cleft*alpha*tend))^(2/(gamma-1));
    uexact(i) = uleft + (2/(gamma+1))*(xx(i)-pos1)/tend;
    cexact(i) = sqrt(gamma*pexact(i)/rhoexact(i));
    machexact(i) = uexact(i)/cexact(i);
  elseif xx(i) <= conpos
    pexact(i) = p3;    	      rhoexact(i) = rho2;
    uexact(i) = u2+uright;    cexact(i)   = sqrt(gamma*pexact(i)/rhoexact(i));
    machexact(i) = uexact(i)/cexact(i);
  elseif xx(i) <= spos
    pexact(i) = p3;    rhoexact(i) = rho3;    uexact(i) = u2+uright; 
    cexact(i) = sqrt(gamma*pexact(i)/rhoexact(i));
    machexact(i) = uexact(i)/cexact(i);
  else
    pexact(i) = pright;    rhoexact(i) = rhoright;
    uexact(i) = uright;    cexact(i)   = sqrt(gamma*pexact(i)/rhoexact(i));
    machexact(i) = uexact(i)/cexact(i);
  end
end
entroexact = log(pexact./rhoexact.^gamma);
energexact = gammab*(pexact./rhoexact);

% ---------END OF CODE----------- %