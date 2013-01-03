% Numerical solution for wave propagation

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.3 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Equation (9.34) is solved with the kappa scheme
% dy/dt + c dy/dx = 0
% Known exact solution f(x-ct)
% Boundary conditions using upwind schemes or extrapolation 
% Uniform vertex-centred grid
% This program generates Figs. 9.6 and 9.7 in the book

% Functions called: exact_solution

		% ...............Input.........................................
kappa = 0;	% Parameter of kappa-scheme 
c = 1;		% Velocity  
sigma = 0.7;	% CFL number	 
J = 31;		% Number of grid cells
bc = 2;		% Governs choice of numerical boundary conditions; see page 350
		% Enter 1 for boundary conditions using upwind schemes
		%    or 2 for extrapolation using  equation (9.36)
T = 0.4;	% Final time
		% ..............End of input...................................

h = 1.0/(J-1);		% Size of grid cells 
dt = sigma*h/c;		% Time step
n = floor(T/dt);	% Number of time steps

% 	Definition of grid numbering 
		
%     x=0    			  x=1
% grid |------|------|------|------|
%      1      2      3             J

x = [0:J-1]'*h;		% Location of grid nodes 

sig2 = sigma^2; sig4 = sigma/4; 

ynew = zeros(J,1); yright = zeros(n+1,1);  
t = 0;	yright(1) = exact_solution(t,1,c);  yold = exact_solution(t,x,c);
for i = 1:n,  t = t + dt;
  a = 2*yold(1) - yold(2);	% Extrapolated inflow value
  ynew(1) = exact_solution(t,0,c);
  if bc == 1
    ynew(2) =  yold(2) - sigma*(yold(2) - yold(1));
  else
    ynew(2) =  yold(2) +sig4*(-1+kappa + (1-3*kappa)*sigma +2*kappa*sig2)*a +...
    sig4*(5-3*kappa+(9*kappa-1)*sigma-6*kappa*sig2)*yold(1)...
    +sig4*(-3+3*kappa -(1+9*kappa)*sigma+6*kappa*sig2)*yold(2)+...
    sig4*(-1-kappa+(1+3*kappa)*sigma-2*kappa*sig2)*yold(3);
  end
  for j = 3:J-1
    ynew(j) = yold(j) +sig4*(-1+kappa + (1-3*kappa)*sigma +2*kappa*sig2)*...
    yold(j-2)+sig4*(5-3*kappa+(9*kappa-1)*sigma-6*kappa*sig2)*yold(j-1)...
    +sig4*(-3+3*kappa -(1+9*kappa)*sigma+6*kappa*sig2)*yold(j)+...
    sig4*(-1-kappa+(1+3*kappa)*sigma-2*kappa*sig2)*yold(j+1);
  end
  if bc == 1		% Fully upwind scheme at outflow
    ynew(J) = (1-1.5*sigma+0.5*sig2)*yold(J) + sigma*(2-sigma)*yold(J-1)...
    +sigma*0.5*(sigma-1)*yold(J-2);      
  else			% Extrapolation at outflow
      ynew(J) = yold(J) +sig4*(-1+kappa + (1-3*kappa)*sigma +2*kappa*sig2)*...
    yold(J-2)+sig4*(5-3*kappa+(9*kappa-1)*sigma-6*kappa*sig2)*yold(J-1)...
    +sig4*(-3+3*kappa -(1+9*kappa)*sigma+6*kappa*sig2)*yold(J)+... 
      sig4*(-1-kappa+(1+3*kappa)*sigma-2*kappa*sig2)*(2*yold(J) - yold(J-1));
  end
  yright(i+1) = ynew(J);  yold = ynew;
end

yexact =  exact_solution(t,x,c);	error = ynew - yexact;
norm1 = norm(error,1)/J			% Compute error norms
norm2 = norm(error,2)/sqrt(J)
norminf = norm(error,inf)

figure(1), clf, hold on
xx = 0:0.005:1;	plot (xx,exact_solution(t,xx,c),'-'); plot (x,ynew,'o');
s1 = ['kappa=',num2str(kappa),' bc =', num2str(bc),' sigma=',...
	num2str(sigma),' cells=',num2str(J-1),' t=',num2str(t)];
title(s1,'fontsize',18);
