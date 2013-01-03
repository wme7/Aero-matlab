% Multistage omega scheme.

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 5.10 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program generates Figs. 5.23, 5.24 in the book

% Solution of (5.89)
% dy/dt + u dy/dx - D d^2u/dx^2 = (hom-1)*beta^2*D*cos(beta*(x-ut))
%	x = 0: Dirichlet,  x = 1: Neumann or absorbing boundary conditions
% Central or upwind differences
% Defect correction optional
% Uniform cell-centered grid

% Functions called: exact_solution, fL, fR

		% ...............Input.........................................
u = 1.1;	% Velocity 
D = 0.0005;	% Diffusion coefficient	
dt = 1/10;	% Time step
n = 100;	% Number of time steps
J = 30;		% Number of grid cells
central = 1;	% Enter 1 for central scheme or something else for upwind scheme
bc = 2;		% Enter 1 for inhomogeneous Neumann or 2 for homogeneous Neumann
hom = 1;	% Enter 1 for   homogeneous right-hand side 
		%    or 2 for inhomogeneous right-hand side
defcor = 1;	% Enter 1 for defect correction or else something else
nalpha = 3;	% Parameter in exact solution: alpha = nalpha*pi
nbeta = 2;	% Parameter in exact solution: beta  = nbeta*pi
		%...............End of input...............................

omega = 1 - sqrt(0.5);			% Parameters in multistage omega scheme
palpha = (1 - 2*omega)/(1 - omega);
palphap = 1 - palpha;
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

	% Preallocation of J*J tridiagonal matrix A
Ac = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);  % Central
Au = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);  % Upwind

% Central scheme
Ac(1,1) = sigma/2 + 3*d;  Ac(1,2) = sigma/2 - d;
for j = 2:J-1,
  Ac(j,j-1) = -sigma/2 - d;    Ac(j,j+1) = sigma/2 - d;
  Ac(j,j)   = - Ac(j,j-1) - Ac(j,j+1);
end
j = J; Ac(j,j-1) = -sigma/2 - d; Ac(j,j) = sigma/2 + d;

% Upwind scheme
Au(1,1) = sigma + 3*d;  Au(1,2) =  - d;
for j = 2:J-1,
  Au(j,j-1) = -sigma - d;  Au(j,j+1) = - d;
  Au(j,j)   = - Au(j,j-1) - Au(j,j+1);
end
j = J; Au(j,j-1) = -sigma - d;  Au(j,j) = sigma + d;

	% System to be solved is E*ynew = C*yold + rhs
z = zeros(J,1); e = ones(J,1);   ID = spdiags([z e z],[-1 0 1],J,J);
Ec1 = ID + palpha*omega*Ac;	Cc1 = ID - palphap*omega*Ac;
Ec2 = ID + palphap*omegap*Ac;	Cc2 = ID - palpha*omegap*Ac;
Eu1 = ID + palpha*omega*Au;	Cu1 = ID - palphap*omega*Au;
Eu2 = ID + palphap*omegap*Au;	Cu2 = ID - palpha*omegap*Au;

if defcor ~= 1		% No defect correction
  if central == 1		% Central scheme
    E1 = Ec1;    C1 = Cc1;    E2 = Ec2;    C2 = Cc2;
  else				% Upwind scheme
    E1 = Eu1;    C1 = Cu1;    E2 = Eu2;    C2 = Cu2;
  end
end

	% Preallocations
yold = zeros(J,1);		% Numerical solution at previous time
ynew = zeros(J,1);		% Numerical solution at new time

x = [1:J]'*h - h/2;		% Coordinates of cell centers

%for j = 1:J,  x(j)=(j-0.5)*h; end

t = 0; yold = exact_solution(t,x,D,nalpha,nbeta,u,hom-1);
rhs = zeros(J,1);		% Right-hand side
beta = nbeta*pi;

if defcor ~= 1		% No defect correction
for i = 1:n

% Stage 1
  rhs = (hom-1)*beta*beta*D*dt*omega*cos(beta*(x-u*t));
  rhs(1) = rhs(1) +  palphap*omega*(sigma + 2*d)*fL(t,D,nalpha,nbeta,u,hom-1)...
   + palpha*omega*(sigma + 2*d)*fL(t+omega*dt,D,nalpha,nbeta,u,hom-1);
  if central == 1		% Central scheme
    rhs(J) = rhs(J) -...
     palphap*omega*(sigma/2 - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
     - palpha*omega*(sigma/2 - d)*fR(t+omega*dt,D,nalpha,nbeta,u,bc,hom-1);
  else				% Upwind scheme
    rhs(J) = rhs(J) - palphap*omega*( - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
     - palpha*omega*( - d)*fR(t+ omega*dt,D,nalpha,nbeta,u,bc,hom-1);
  end
  rhs = rhs + C1*yold;  ynew = E1\rhs;

% Stage 2
  yold = ynew;
  rhs = (hom-1)*beta*beta*D*dt*omegap*cos(beta*(x-u*(t + dt - omega*dt)));
  rhs(1) = rhs(1) + ...
     palpha*omegap*(sigma + 2*d)*fL(t + omega*dt,D,nalpha,nbeta,u,hom-1)...
   + palphap*omegap*(sigma + 2*d)*fL(t + dt - omega*dt,D,nalpha,nbeta,u,hom-1);
  if central == 1		% Central scheme
    rhs(J) = rhs(J) -...
     palpha*omegap*(sigma/2 - d)*h*fR(t + omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
     - palphap*omegap*(sigma/2 - d)*fR(t+dt-omega*dt,D,nalpha,nbeta,u,bc,hom-1);
  else				% Upwind scheme
    rhs(J) = rhs(J) -...
      palpha*omegap*( - d)*h*fR(t + omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
     - palphap*omegap*( - d)*fR(t + dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1);
  end
  rhs = rhs + C2*yold;  ynew = E2\rhs;

% Stage 3
  yold = ynew;
  rhs = (hom-1)*beta*beta*D*dt*omega*cos(beta*(x - u*(t + dt)));
  rhs(1) = rhs(1) +...
    palphap*omega*(sigma + 2*d)*fL(t+ dt- omega*dt,D,nalpha,nbeta,u,hom-1)...
   + palpha*omega*(sigma + 2*d)*fL(t+dt,D,nalpha,nbeta,u,hom-1);
  if central == 1		% Central scheme
    rhs(J) = rhs(J) -...
    palphap*omega*(sigma/2-d)*h*fR(t+ dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
     - palpha*omega*(sigma/2 - d)*fR(t+dt,D,nalpha,nbeta,u,bc,hom-1);
  else				% Upwind scheme
    rhs(J) = rhs(J) -...
   palphap*omega*( - d)*h*fR(t+ dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
     - palpha*omega*( - d)*fR(t+ dt,D,nalpha,nbeta,u,bc,hom-1);
  end
  rhs = rhs + C1*yold;  ynew = E1\rhs;
  t = t + dt;  yold = ynew;
end
else		% One step of defect correction
  dy = zeros(J,1);		% Correction
  drhs = zeros(J,1);		% Correction of right-hand side
  
  for i = 1:n
		% Stage 1: upwind step
    rhs = (hom-1)*beta*beta*D*dt*omega*cos(beta*(x-u*t));
    rhs(1) = rhs(1) +  palphap*omega*(sigma + 2*d)*fL(t,D,nalpha,nbeta,u,hom-1)...
     + palpha*omega*(sigma + 2*d)*fL(t+omega*dt,D,nalpha,nbeta,u,hom-1);

    drhs = rhs;		% Preparation of rhs for correction step
    drhs(J) = drhs(J) -...
     palphap*omega*(sigma/2 - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
     - palpha*omega*(sigma/2 - d)*fR(t+omega*dt,D,nalpha,nbeta,u,bc,hom-1);
    rhs(J) = rhs(J) - palphap*omega*( - d)*h*fR(t,D,nalpha,nbeta,u,bc,hom-1)...
       - palpha*omega*( - d)*fR(t+ omega*dt,D,nalpha,nbeta,u,bc,hom-1);
    rhs = rhs + Cu1*yold;    ynew = Eu1\rhs;
    
		% Stage 1: correction step
    drhs = drhs - Ec1*ynew + Cc1*yold;    dy = Eu1\drhs;    ynew = ynew + dy;

		% Stage 2: upwind step
    yold = ynew;
    rhs = (hom-1)*beta*beta*D*dt*omegap*cos(beta*(x-u*(t + dt - omega*dt)));
    rhs(1) = rhs(1) + ...
      palpha*omegap*(sigma + 2*d)*fL(t + omega*dt,D,nalpha,nbeta,u,hom-1)...
    + palphap*omegap*(sigma + 2*d)*fL(t + dt - omega*dt,D,nalpha,nbeta,u,hom-1);

    drhs = rhs;		% Preparation of rhs for correction step
    drhs(J) = drhs(J) -...
       palpha*omegap*(sigma/2 - d)*h*fR(t + omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
       - palphap*omegap*(sigma/2 - d)*fR(t+dt-omega*dt,D,nalpha,nbeta,u,bc,hom-1);
    rhs(J) = rhs(J) -...
       palpha*omegap*( - d)*h*fR(t + omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
       - palphap*omegap*( - d)*fR(t + dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1);
    rhs = rhs + Cu2*yold;    ynew = Eu2\rhs;
    
		% Stage 2: correction step
    drhs = drhs - Ec2*ynew + Cc2*yold;    dy = Eu2\drhs;    ynew = ynew + dy;

		% Stage 3: upwind step
    yold = ynew;
    rhs = (hom-1)*beta*beta*D*dt*omega*cos(beta*(x - u*(t + dt)));
    rhs(1) = rhs(1) +...
    palphap*omega*(sigma + 2*d)*fL(t+ dt- omega*dt,D,nalpha,nbeta,u,hom-1)...
    + palpha*omega*(sigma + 2*d)*fL(t+dt,D,nalpha,nbeta,u,hom-1);

    drhs = rhs;		% Preparation of rhs for correction step
    drhs(J) =drhs(J) -...
      palphap*omega*(sigma/2-d)*h*fR(t+ dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
       - palpha*omega*(sigma/2 - d)*fR(t+dt,D,nalpha,nbeta,u,bc,hom-1);
      rhs(J) = rhs(J) -...
      palphap*omega*( - d)*h*fR(t+ dt- omega*dt,D,nalpha,nbeta,u,bc,hom-1)...
       - palpha*omega*( - d)*fR(t+ dt,D,nalpha,nbeta,u,bc,hom-1);
    rhs = rhs + Cu1*yold;
    ynew = Eu1\rhs;
    
		% Stage 3: correction step
    drhs = drhs - Ec1*ynew + Cc1*yold;    dy = Eu1\drhs;    ynew = ynew + dy;
    t = t + dt;    yold = ynew;
  end
end

solex =  exact_solution(t,x,D,nalpha,nbeta,u,hom-1);

error = ynew - solex;		% Compute error norms
norm1 = norm(error,1)/J
norm2 = norm(error,2)/sqrt(J)
norminf = norm(error,inf)

figure(1), clf, hold on
xx = 0:0.02:1; plot (xx,exact_solution(t,xx,D,nalpha,nbeta,u,hom-1),'-');
axis([0 1 -1 1]); plot (x,ynew,'o');
if central == 1,  s1 = ['Central fractional step omega-scheme'];
else,             s1 = ['Upwind fractional step omega-scheme'];
end
if defcor == 1
  s1 = ['Fractional step omega-scheme with defect correction'];
end
if bc == 1,  s1 = [s1, ', inhomogeneous Neumann'];
else         s1 = [s1, ', homogeneous Neumann'];
end
title(s1,'fontsize',12);
if hom ==1,  s2 = ['Homogeneous differential equation']
else,        s2 = ['Inhomogeneous differential equation']
end
s3 = ['omega=',num2str(omega),'  D=',num2str(D),' u=',num2str(u),...
 ' h=',num2str(h),'  dt=',num2str(dt),'  t=',num2str(T),'  n=',num2str(n)]
s4 = ['meshpeclet=',num2str(meshpeclet),' Courant=',num2str(sigma),...
 '  d=',num2str(2*d) ]
