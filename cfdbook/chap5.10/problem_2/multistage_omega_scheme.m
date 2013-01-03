% Multistage omega scheme.

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 5.10 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Fig. 5.22 in the book

% Solution of (5.89)
% dy/dt + u dy/dx - D d^2u/dx^2 = 0
% Initial and boundary conditions are such that with D = 0 we would have a
% propagating block function
%	x = 0: Dirichlet,  x = 1: homogeneous Neumann boundary condition
% Central or upwind differences
% Uniform cell-centered grid

% Functions called:  fL, fstart
% fL has to be adapted, depending on the velocity u being zero or not

		% ...............Input.........................................
u = 0.0;	% Velocity 
D = 1/50;	% Diffusion coefficient
dt = 1/5;	% Time step
n = 1;		% Number of time steps
J = 50;		% Number of grid cells
central = 1;	% Enter 1 for central scheme or something else for upwind scheme
		%...............End of input...............................

omega = 1 - sqrt(0.5);			% Parameters in multistage omega scheme
alpha = (1 - 2*omega)/(1 - omega);
alphap = 1 - alpha;
omegap = 1 - 2*omega;

T = n*dt;	% Final time

		% Definition of indexing  
%     x=0    
% grid |---o---|---o---|---o---|---o---|
%      1   1   2   2   3           J

h = 1.0/J;			% Cell size
d = D*dt/(h*h);			% Diffusion number
sigma = u*dt/h;			% Courant number
meshpeclet = sigma/d;		% Mesh Peclet number

	% Declaration of J*J tridiagonal matrix A
A = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);  

if central == 1  	% Central scheme
  A(1,1) = sigma/2 + 3*d;       	A(1,2)   = sigma/2 - d;
  for j = 2:J-1,
    A(j,j-1) = -sigma/2 - d;    	A(j,j+1) = sigma/2 - d;
    A(j,j)   = - A(j,j-1) - A(j,j+1);
  end
  j = J;  A(j,j-1) = -sigma/2 - d;  	A(j,j)   = sigma/2 + d;
else  			% Upwind scheme
  A(1,1) = sigma + 3*d;  A(1,2) =  - d;
  for j = 2:J-1,
    A(j,j-1) = -sigma - d;	        A(j,j+1) = - d;
    A(j,j)   = - A(j,j-1) - A(j,j+1);
  end
  j = J;  A(j,j-1) = -sigma - d;  	A(j,j)   = sigma + d;
end

	% System to be solved is E*ynew = C*yold + rhs
z = zeros(J,1); e = ones(J,1); ID = spdiags([z e z],[-1 0 1],J,J);
E1 = ID + alpha*omega*A;	C1 = ID - alphap*omega*A;
E2 = ID + alphap*omegap*A;	C2 = ID - alpha*omegap*A;

	%Preallocations
yold = zeros(J,1);		% Numerical solution at previous time
ynew = zeros(J,1);		% Numerical solution at new time

x = [1:J]'*h - h/2;		% Coordinates of cell centers

for j = 1:J,  yold(j) = fstart(x(j)); end

	% Plot of initial condition
figure(1), clf, hold on, plot (x,yold); axis([0 1 -1.5 1.5]);

t = 0;
for i = 1:n

% Stage 1
  rhs = zeros(J,1);
  rhs(1) = rhs(1) +  alphap*omega*(sigma + 2*d)*fL(u*t)...
   + alpha*omega*(sigma + 2*d)*fL(u*(t+omega*dt));
  rhs = rhs + C1*yold;  ynew = E1\rhs;  
  
% Stage 2
  yold = ynew; 	rhs = zeros(J,1);
  rhs(1) = rhs(1) +  alpha*omegap*(sigma + 2*d)*fL(u*(t+omega*dt))...
   + alphap*omegap*(sigma + 2*d)*fL(u*(t+ dt - omega*dt));
  rhs = rhs + C2*yold;  ynew = E2\rhs;

% Stage 3
  yold = ynew;  rhs = zeros(J,1);
  rhs(1) = rhs(1) +  alphap*omega*(sigma + 2*d)*fL(u*(t+ dt - omega*dt))...
   + alpha*omega*(sigma + 2*d)*fL(u*(t+dt));
  rhs = rhs + C1*yold;  ynew = E1\rhs;
  
  t = t + dt;  yold = ynew;
end

plot (x,ynew,'o');
if central == 1,  	s1 = ['Central fractional step omega-scheme'];
else,  			s1 = ['Upwind fractional step omega-scheme'];
end
title(s1,'fontsize',16);
s3 = ['omega=',num2str(omega),' alpha =', num2str(alpha),' D=',num2str(D),...
  ' u=',num2str(u),' h=',num2str(h),'  dt=',num2str(dt)]
s4 = ['  t=',num2str(t),' meshpeclet=',num2str(meshpeclet),' Courant=',...
  num2str(sigma),'  d=',num2str(2*d),'  n=',num2str(n) ]
