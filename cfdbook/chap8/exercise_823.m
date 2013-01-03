% Exercise 8.2.3

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Numerical solution of 1-D shallow water equations (8.11), (8.12):
%	dh/dt + h*.du/dx = 0
%	du/dt + g.dh/dx = 0

% Theory is given in Section 8.2 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
%	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program solves Exercise 8.2.3 
% Uniform staggered grid (Fig. 8.2)
% x contains h-points; y contains u-points; xy contains all points

% Functions called: hinit ,uinit

clear all 

%****************************INPUT********************************************

hstar = 0.2; g = 1;	% Coefficients in shallow water equations
scheme = 4;		% 1: Implicit scheme (5.76); unconditionally stable
			% 2: Hansen scheme (5.59,5.60); stable for sigma <= 1  
			% 3: Leapfrog scheme (5.48); stable for sigma <= 0.5
			% 4: Sielecki scheme (5.75); stable for sigma <= 1 
geval = 3; 		% 1: Piecewise linear initial depth
			% 2: Sinusoidal initial depth
			% 3: Hyperbolic tangent initial depth 
sigma = 0.4;		% Courant number c*dt/dx; c = sqrt(g*hstar)
m = 20;			% Number of cells
tend = 0.6;		% Final time

%****************************END OF INPUT************************************

%		Definition of grid                              
%     	x=0       m cells       x=1
% 	grid |---o---|---o---|---o---|
% 	     1   1   2   2   3   3   4
% 	     h   u   h   u   h   u   h

dx = 1/m; x = dx*[0:1:m];		% x-coordinates of h-nodes
y = dx*[1:1:m] - 0.5*dx;		% x-coordinates of u-nodes

	% Preallocations
solold = zeros(m+m+1,1); solnew = solold;  % Solution: sol = (h1,u1,h2,u2,....)
h = zeros(m+1,1); hold = h; 		% h = (h1,h2,...)
u = zeros(m,1); uold = u; 		% u = (u1,u2,...)

c = sqrt(g*hstar);			% c is the wave velocity
dt = sigma*dx/c; % Timestep follows from the definition of the Courant number
ntime = floor(tend/dt); 		% Number of time steps

for j = 1:m+1,  hold(j) = hstar*hinit(x(j),geval); end	% Initial h 
for j = 1:m,    uold(j) = c*uinit(y(j));	   end	% Initial u
 
if scheme == 1	% Implicit scheme (5.76)

  for j = 1:2:m+m+1
    z = (j-1)*dx/2; solold(j) = hstar*hinit(z,geval);	% Initial h
  end
  for j = 2:2:m+m
    z = (j-1)*dx/2; solold(j) = c*uinit(z);		% Initial u
  end
	% Preallocation of J*J tridiagonal matrix A
  J = m + m + 1;
  A = spdiags([ ones(J,1)  ones(J,1)  ones(J,1)], [-1 0 1]',J,J);
	% Matrix generation
  htx = hstar*dt/(2*dx); gtx = g*dt/(2*dx);
  A(1,1) = 1;
  for j = 2:2:m+m
    A(j,j-1) = -gtx; A(j,j) = 1; A(j,j+1) = gtx;
  end
  for j = 3:2:m+m-1
    A(j,j-1) = -htx; A(j,j) = 1; A(j,j+1) = htx;
  end
  j = m+m+1; A(j,j) = 1;

  f = zeros(m+m+1,1);		% Preallocation of right-hand side
  t=0;				% Time-stepping
  for i=1:ntime,   t = t+dt;
    f(1) = 0.5*hstar*(hinit(-c*t,geval)+uinit(-c*t)+hinit(c*t,geval)...
    -uinit(c*t));
    for j = 2:2:m+m
      f(j) = solold(j) - gtx*(solold(j+1)-solold(j-1));
    end
    for j = 3:2:m+m-1
      f(j) =  solold(j) - htx*(solold(j+1)-solold(j-1));
    end
    f(m+m+1)=0.5*hstar*(hinit(1-c*t,geval)+uinit(1-c*t)+hinit(1+c*t,geval)...
    -uinit(1+c*t));
    solnew = A\f; solold = solnew;
  end

  for j = 1:m+1,  h(j) = solold(2*j -1); end	% Extraction of h and u
  for j = 1:m,    u(j) = solold(2*j);    end
end;			% End of implicit scheme

if scheme == 2 	 	% Hansen scheme (8.34), (8.35)
%  Initial solution has been computed above.
%  Because Hansen scheme is staggered in time we use exact solution for initial
%  velocity at t = dt/2
  for j = 1:m
    z = y(j); uold(j) = 0.5*c*(hinit(z-c*dt/2,geval)+uinit(z-c*dt/2)-...
              hinit(z+c*dt/2,geval) + uinit(z+c*dt/2));
  end
  htx = hstar*dt/dx; gtx = g*dt/dx;
  for i = 1:ntime, t = i*dt;
	% We prescribe the exact solution at the boundaries
    h(1) = 0.5*hstar*(hinit(x(1)-c*t,geval)+uinit(x(1)-c*t)+...
    hinit(x(1)+c*t,geval) - uinit(x(1)+c*t));
    h(m+1) = 0.5*hstar*(hinit(x(m+1)-c*t,geval)+uinit(x(m+1)-c*t)+...
    hinit(x(m+1)+c*t,geval) - uinit(x(m+1)+c*t));
    for j = 2:m
      h(j) = hold(j) - htx*(uold(j) - uold(j-1));
    end
    for j = 1:m
      u(j) = uold(j) - gtx*(h(j+1) - h(j));
    end      
    uold = u; hold = h;         
  end
end		% End of Hansen scheme

if scheme == 3	% Leapfrog scheme (8.27). This scheme requires extra initial
		% conditions that are taken from the exact solution. 
  hveryold = hold; 		% Solution has to be stored at two previous
  uveryold = uold; 		% time levels
  t = dt; htx = hstar*dt/dx; gtx = g*dt/dx;
  for j = 1:m+1,  z = (j-1)*dx;
    hold(j) = 0.5*hstar*(hinit(z-c*t,geval)+uinit(z-c*t)+...
    hinit(z+c*t,geval) - uinit(z+c*t));
  end
  for j = 1:m,  z = (j-1/2)*dx;
    uold(j) = 0.5*c*(hinit(z-c*t,geval)+uinit(z-c*t)-...
    hinit(z+c*t,geval) + uinit(z+c*t));
  end
  for i = 2:ntime
    t = i*dt;		% We prescribe the exact solution at the boundaries
    h(1) = 0.5*hstar*(hinit(x(1)-c*t,geval)+uinit(x(1)-c*t)+...
    hinit(x(1)+c*t,geval) - uinit(x(1)+c*t));
    h(m+1) = 0.5*hstar*(hinit(x(m+1)-c*t,geval)+uinit(x(m+1)-c*t)...
    + hinit(x(m+1)+c*t,geval) - uinit(x(m+1)+c*t));
    for j = 2:m
      h(j) = hveryold(j) - 2*htx*(uold(j) - uold(j-1));
    end
    for j = 1:m
      u(j) = uveryold(j) - 2*gtx*(hold(j+1) - hold(j));
    end      
    uveryold = uold; hveryold = hold; uold = u; hold = h;         
  end  
end 			% End of leapfrog scheme

if scheme == 4 	 	% Sielecki scheme.
			% Initial solution has been computed above
  htx = hstar*dt/dx; gtx = g*dt/dx;
  for i = 1:ntime
    t = i*dt;	% We prescribe the exact solution at the boundaries
    h(1) = 0.5*hstar*(hinit(x(1)-c*t,geval)+uinit(x(1)-c*t)+...
    hinit(x(1)+c*t,geval) - uinit(x(1)+c*t));
    h(m+1) = 0.5*hstar*(hinit(x(m+1)-c*t,geval)+uinit(x(m+1)-c*t)+...
    hinit(x(m+1)+c*t,geval) - uinit(x(m+1)+c*t));
    for j = 2:m
      h(j) = hold(j) - htx*(uold(j) - uold(j-1));
    end
    for j = 1:m
      u(j) = uold(j) - gtx*(h(j+1) - h(j));
    end      
    uold = u; hold = h;         
  end
end			% End of Sielecki scheme

	% Exact solution
zz = 0:0.01:1; hexact = zeros(length(zz),1); uexact = hexact;
for j = 1:length(zz),
  hexact(j) = 0.5*hstar*(hinit(zz(j)-c*t,geval)+uinit(zz(j)-c*t)+...
  hinit(zz(j)+c*t,geval) - uinit(zz(j)+c*t));

  if scheme == 2	% Hansen scheme is staggered in time
    cctt = c*(t+dt/2);
    uexact(j) = 0.5*c*(hinit(zz(j)-cctt,geval)+uinit(zz(j)-cctt)-...
    hinit(zz(j)+cctt,geval) + uinit(zz(j)+cctt));
  else
    uexact(j) = 0.5*c*(hinit(zz(j)-c*t,geval)+uinit(zz(j)-c*t)-...
    hinit(zz(j)+c*t,geval) + uinit(zz(j)+c*t));
  end
end

figure(1), clf, subplot(211), hold on
plot(zz,hexact)
if     scheme == 1, text = ' Implicit';
elseif scheme == 2, text = ' Hansen';
elseif scheme == 3, text = ' Leapfrog';
else, 		    text = ' Sielecki';
end
title(['Water Level',text,' Courant = ',num2str(sigma),', ',num2str(m),...
' cells, time: ', num2str(t)])
plot(x,h,'*')
subplot(212), hold on
plot(zz,uexact), plot(y,u,'*'), title(['Velocity'])
