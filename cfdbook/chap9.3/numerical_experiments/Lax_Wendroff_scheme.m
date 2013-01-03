% Numerical solution for wave propagation

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 9.3 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Equation (9.34) is solved with the Lax-Wendroff scheme
% dy/dt + c dy/dx = 0
% Known exact solution f(x-ct)
% Boundary conditions using upwind schemes or extrapolation 
% Uniform vertex-centred grid
% This program generates Figs. 9.5 in the book

% Functions called: exact_solution

		% ...............Input.........................................
c = 1;		% Velocity
sigma = 0.7;	% CFL number 
J = 31;		% Number of grid cells
bc = 1;		% Governs choice of numerical boundary conditions; see page 350
		% Enter 1 for boundary conditions using upwind schemes
		%    or 2 for extrapolation using  equation (9.36)
T = 0.4;	% Final time
		% ..............End of input...................................

h = 1.0/(J-1);		% Size of grid cells  
dt = sigma*h/c;		% Time step
n = floor(T/dt);	% Number of time steps

% 		  Definition of grid numbering
%	     x=0    			 	  x=1
%	 grid |------|------|---.........---|------|
%	      1      2      3             	   J

x = [0:J-1]'*h;		% Location of grid nodes 

sig2 = sigma^2;
t = 0;
yold = exact_solution(t,x,c); 	ynew = zeros(J,1);
yright = zeros(n+1,1);  	yright(1) = exact_solution(t,1,c);

for i = 1:n,  t = t + dt;  ynew(1) = exact_solution(t,0,c);
  for j = 2:J-1
    ynew(j) = 0.5*(sig2+sigma)*yold(j-1)+(1-sig2)*yold(j)+...
     0.5*(sig2-sigma)*yold(j+1);
  end
  if bc == 1,		% First order upwind at outflow
     ynew(J) = yold(J) + sigma*(yold(J-1) - yold(J));
  else			% Extrapolation at outflow
     ynew(J) = 0.5*(sig2+sigma)*yold(J-1)+(1-sig2)*yold(J)+...
     0.5*(sig2-sigma)*(2*yold(J) - yold(J-1));
  end
  yright(i+1) = ynew(J);  yold = ynew;
end

yexact =  exact_solution(t,x,c);	error = ynew - yexact;
norm1 = norm(error,1)/J			% Compute error norms
norm2 = norm(error,2)/sqrt(J)
norminf = norm(error,inf)

figure(1), clf, hold on
xx = 0:0.005:1;	plot (xx,exact_solution(t,xx,c),'-'); plot (x,ynew,'o');
s1 = ['Lax-Wendroff',' bc =', num2str(bc),' sigma=',...
num2str(sigma),' cells=',num2str(J-1),' t=',num2str(t)];
title(s1,'fontsize',18);





