% Numerical experiments in Section 5.9 of

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% This program generates the following figures in the book:
% left part of Figs. 5.17 -- 5.21

% Solution of (5.89)
% dy/dt + u dy/dx - D d^2u/dx^2 = 0
%   x = 0: Dirichlet, x = 1: homogeneous Neumann or absorbing boundary condition
% Step function as initial solution
% Central or upwind differences
% Uniform cell-centered grid

% Functions called: fL, fstart
% See body of fL for modification in case of Figs 5.17 and 5.18

		% ...............Input.........................................
u = 0.9;	% Velocity 
D = 0.002;	% Diffusion coefficient
dt = 1/50;	% Time step
n = 20;	% Number of time steps
J = 50;		% Number of grid cells
omega = 0.5;	% Parameter of omega time stepping method
central = 2;	% Enter 1 for central scheme or something else for upwind scheme
bc = 1;		% Enter 2 for transparent boundary condition or something else  
		% for homogeneous Neumann boundary condition
		%...............End of input...............................

T = n*dt;	% Final time

		% Definition of indexing 
%     x=0    
% grid |---o---|---o---|---o---|---o---|
%      1   1   2   2   3           J

h = 1.0/J;			% Cell size
d = D*dt/(h*h);			% Diffusion number
sigma = u*dt/h;			% Courant number
meshpeclet = sigma/d;		% Mesh Peclet number

	% Preallocation of J*J tridiagonal matrix A
A = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);

if central == 1	  		% Central scheme
  A(1,1) = sigma/2 + 3*d;  A(1,2) = sigma/2 - d;
  for j = 2:J-1,
    A(j,j-1) = - sigma/2 - d;    A(j,j+1) =  sigma/2 - d;
    A(j,j)   = - A(j,j-1) - A(j,j+1);
  end
  j = J;
  if bc == 2		  	% Transparent boundary condition
    A(j,j-1) = -sigma;    A(j,j)   =  sigma;
  else				% Homogeneous Neumann boundary condition
    A(j,j-1) = -sigma/2 - d;    A(j,j)   =  sigma/2 + d;
  end
else				% Upwind scheme
  A(1,1) = sigma + 3*d;  A(1,2) =  - d;
  for j = 2:J-1,
    A(j,j-1) = - sigma - d;    A(j,j+1) = - d;
    A(j,j)   = - A(j,j-1) - A(j,j+1);
  end
  j = J;
  if bc == 2		  	% Transparent boundary condition
    A(j,j-1) = -sigma;    A(j,j)   =  sigma;
  else				% Homogeneous Neumann boundary condition
    A(j,j-1) = -sigma - d;    A(j,j)   =  sigma + d;
  end
end

E = eye(size(A)) + omega*A;	% System to be solved is E*ynew = C*yold + rhs 
C = eye(size(A)) + (omega-1)*A;		
yold = zeros(J,1);		% Numerical solution at previous time
ynew = zeros(J,1);		% Numerical solution at new time
x = zeros(J,1);			% Coordinates of cell centers
for j = 1:J
  x(j) = (j-0.5)*h;  yold(j) = fstart(x(j));
end

figure(1), clf, hold on, plot (x,yold);	% Plot initial condition
axis([0 1 -1.5 1.5]);

t = 0;
for i = 1:n
  rhs = zeros(J,1);
  rhs(1) = rhs(1) + (1 - omega)*(sigma + 2*d)*fL(u*t)...
   + omega*(sigma + 2*d)*fL(u*(t+dt));
  rhs = rhs + C*yold;  ynew = E\rhs;  t = t + dt;  yold = ynew;
end
plot (x,ynew,'o');
if central == 1,  s1 = ['Central omega scheme'];
else,             s1 = ['Upwind omega scheme'];
end
if bc == 1,  	  s1 = [s1, ', homogeneous Neumann'];
else,       	  s1 = [s1, ', absorbing'];
end
title(s1,'fontsize',18)
s3 = ['omega=',num2str(omega),'  D=',num2str(D),' u=',num2str(u),...
 ' h=',num2str(h),'  dt=',num2str(dt),'  t=',num2str(T),'  n=',num2str(n)]
s4 = ['meshpeclet=',num2str(meshpeclet),' Courant=',num2str(sigma),...
 '  d=',num2str(2*d) ]
