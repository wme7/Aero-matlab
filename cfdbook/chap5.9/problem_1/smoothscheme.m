% Numerical experiments in Section 5.9 of

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% This program generates the following figures in the book:
% left part of Figs. 5.11, 5.12, 5.14, 5.15, and Figs. 5.13 and 5.16

% Solution of (5.89)
% dy/dt + u dy/dx - D d^2u/dx^2 = (hom-1)*beta^2*D*cos(beta*(x-ut))
%	x = 0: Dirichlet  x = 1: Neumann or absorbing boundary conditions
% Central or upwind differences
% Uniform cell-centered grid

% Functions called: exact_solution, fL, fR

		% ...............Input.........................................
u = 1.1;	% Velocity 
D = 0.02;	% Diffusion coefficient
dt = 1/30;	% Time step
n = 30;		% Number of time steps	
J = 30;		% Number of grid cells
nalpha = 4;	% Parameter in exact solution: alpha = nalpha*pi
nbeta = 2;	% Parameter in exact solution: beta  = nbeta*pi
omega = 1;	% Parameter of omega time stepping method
central = 1;	% Enter 1 for central scheme or something else for upwind scheme
bc = 1;		% Enter 1 for inhomogeneous Neumann or 2 for homogeneous Neumann
                % or 3 for transparent boundary condition
hom = 2;	% Enter 1 for   homogeneous right-hand side 
		%    or 2 for inhomogeneous right-hand side
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

	% Declaration of J*J tridiagonal matrix A
A = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);

if central == 1	  		% Central scheme
  A(1,1) = sigma/2 + 3*d;  A(1,2) = sigma/2 - d;
  for j = 2:J-1,
    A(j,j-1) = - sigma/2 - d;    A(j,j+1) =  sigma/2 - d;
    A(j,j)   = - A(j,j-1) - A(j,j+1);
  end
  j = J;
  if bc == 3		  	% Transparent boundary condition
    A(j,j-1) = -sigma;    A(j,j)   =  sigma;
  else				% Neumann boundary condition
    A(j,j-1) = -sigma/2 - d;    A(j,j)   =  sigma/2 + d;
  end
else				% Upwind scheme
  A(1,1) = sigma + 3*d;  A(1,2) =  - d;
  for j = 2:J-1,
    A(j,j-1) = -sigma - d;   A(j,j+1) = - d;  A(j,j) = - A(j,j-1) - A(j,j+1);
  end
  j = J;
  if bc == 3		  	% Weakly reflecting boundary condition
    A(j,j-1) = -sigma;    A(j,j)   =  sigma;
  else				% Neumann boundary condition
    A(j,j-1) = -sigma - d;    A(j,j)   =  sigma + d;
  end
end

E = eye(size(A)) + omega*A;	% System to be solved is E*ynew = C*yold + rhs
C = eye(size(A)) + (omega-1)*A;
	
yold = zeros(J,1);		% Numerical solution at previous time
ynew = zeros(J,1);		% Numerical solution at new time
x = zeros(J,1);			% Coordinates of cell centers
for j = 1:J,  x(j)=(j-0.5)*h; end

t = 0; yold = exact_solution(t,x,D,nalpha,nbeta,u,hom-1);
rhs = zeros(J,1);		% Right-hand side
beta = nbeta*pi;

for i = 1:n
   rhs = (hom-1)*beta*beta*D*dt*(omega*cos(beta*(x-u*(t+dt)))...
     + (1 - omega)*cos(beta*(x-u*t)));
  rhs(1) = rhs(1) + (1 - omega)*(sigma + 2*d)*fL(t,D,nalpha,nbeta,u,hom-1)...
   + omega*(sigma + 2*d)*fL(t+dt,D,nalpha,nbeta,u,hom-1);
   
  if central == 1		% Central scheme
    rhs(J) = rhs(J) -...
     (1 - omega)*(sigma/2 - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
     - omega*(sigma/2 - d)*fR(t+dt,D,nalpha,nbeta,u,bc,hom-1);
     
    if bc == 3			% Transparent boundary condition    
      rhs(J) = (0.5*sigma/d)*rhs(J);
    end
    
  else				% Upwind scheme
    rhs(J) = rhs(J) - (1 - omega)*( - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
     - omega*( - d)*fR(t+dt,D,nalpha,nbeta,u,bc,hom-1);
  end
  rhs = rhs + C*yold;  ynew = E\rhs;  t = t + dt;  yold = ynew;
end

solex =  exact_solution(t,x,D,nalpha,nbeta,u,hom-1);

error = ynew - solex;		% Compute error norms
norm1 = norm(error,1)/J
norm2 = norm(error,2)/sqrt(J)
norminf = norm(error,inf)

figure(1), clf, hold on
xx = 0:0.02:1; plot (xx,exact_solution(t,xx,D,nalpha,nbeta,u,hom-1),'-');
axis([0 1 -1.5 1.5]); plot (x,ynew,'o');
if central == 1,  s1 = ['Central omega scheme'];
else,             s1 = ['Upwind omega scheme'];
end
if bc == 1,  s1 = [s1, ', inhomogeneous Neumann']; end
if bc == 2,  s1 = [s1, ', homogeneous Neumann'];   end
if bc == 3,  s1 = [s1, ', transparant'];           end
title(s1,'fontsize',18)
if hom == 1,  s2 = ['Homogeneous differential equation']
else,         s2 = ['Inhomogeneous differential equation']
end
s3 = ['omega=',num2str(omega),'  D=',num2str(D),' u=',num2str(u),...
 ' h=',num2str(h),'  dt=',num2str(dt),'  t=',num2str(T),'  n=',num2str(n)]
s4 = ['meshpeclet=',num2str(meshpeclet),' Courant=',num2str(sigma),...
 '  d=',num2str(2*d) ]
s5 = ['nalpha = ',num2str(nalpha),' nbeta =  ',num2str(nbeta)]






