% Numerical solution of Burgers equation: du/dt + 0.5 du^2/dx = 0           #

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory is given in Section 9.4 of
% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
%	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

clear all 

%*****************************INPUT********************************************
geval = 2;	% = 1:     Initial condition: shock: example 9.4.4
		% not = 1: Initial condition: fan:   example 9.4.5
tend = 0.5;	% Final time
lambda = 0.75;	% lambda = dt/dx
J = 40; 	% Number of nodes. 
%******************************************************************************

%*********************Uniform grid*********************
%						      #
%	o------o------o------o------o------o------o   #
%	1      2                                  J   #
%      x=-1                                      x=1  #
%******************************************************

dx = 2/(J-1); dt = lambda*dx; 		% Stepsizes
ntime = floor(tend/dt); %ntime = 0;     % Number of time steps
y = [1:J]*dx - dx -1;			% Location of grid points

u = y;				% Preallocation of solution
if geval == 1, uleft =   1; uright = 0;     % Initial and boundary conditions 
else,  	       uleft = - 1; uright = 1;
end
for j = 1:J,
  if y(j) < -10*eps,    u(j) = uleft;
  elseif y(j) > 10*eps, u(j) = uright;
  else,		  u(j) = 0.5*(uleft+uright);
  end 
end 
uinit = u;

for n = 1:ntime		%   Nonconservative first order upwind scheme (9.75)
  uold = u;
  j = 1; u(j) = uold(j) - 0.5*lambda*...
       ((uold(j) + abs(uold(j)))*  (uold(j) - uleft) + ...
        (uold(j) - abs(uold(j)))*(uold(j+1) - uold(j)));       
  for j = 2:J-1,  u(j) = uold(j) - 0.5*lambda*...
       ((uold(j) + abs(uold(j)))*  (uold(j) - uold(j-1)) + ...
        (uold(j) - abs(uold(j)))*(uold(j+1) - uold(j)));       
  end
  j = J; u(j) = uold(j) - 0.5*lambda*...
     ((uold(j) + abs(uold(j)))*(uold(j) - uold(j-1)) + ...
     (uold(j) - abs(uold(j)))* (uright - uold(j)));
end              
figure(1); clf, subplot(221), hold on, plot(y,u,'o')
title('Nonconservative scheme','fontsize',14)
 
u = uinit;
for n = 1:ntime		% Courant-Isaacson-Rees scheme (9.81)
  uold = u;
  j = 1;  s1 = sign(uold(j) + uold(j+1)); s2 = sign(uold(j) + uleft);
  u(j) = uold(j) - 0.25*lambda*((1+s1)*uold(j)^2 + ...
	    (1-s1)*uold(j+1)^2 - (1+s2)*uleft^2 - (1-s2)*uold(j)^2);
  for j = 2:J-1
    s1   = sign(uold(j) + uold(j+1)); s2 = sign(uold(j) + uold(j-1));
    u(j) = uold(j) - 0.25*lambda*((1+s1)*uold(j)^2 + ...
       (1-s1)*uold(j+1)^2 - (1+s2)*uold(j-1)^2 - (1-s2)*uold(j)^2);
  end
  j = J; s1 = sign(uold(j) + uright); s2 = sign(uold(j) + uold(j-1));
  u(j) = uold(j) - 0.25*lambda*((1+s1)*uold(j)^2 + (1-s1)*uright^2 -...
      (1+s2)*uold(j-1)^2 - (1-s2)*uold(j)^2);
end              
subplot(222), hold on, plot(y,u,'o')
title('Courant-Isaacson-Rees scheme','fontsize',14)

u = uinit;
for n = 1:ntime		% Lax-Friedrichs scheme (9.80)
  uold = u;
  j = 1;      
  u(j) = 0.5*(uleft + uold(j+1)) - 0.25*lambda*(uold(j+1)^2 - uleft^2);
  for j = 2:J-1
    u(j) = 0.5*(uold(j-1)+uold(j+1)) - 0.25*lambda*(uold(j+1)^2 - uold(j-1)^2);
  end
  j = J;      
  u(j) = 0.5*(uold(j-1) + uright) -  0.25*lambda*(uright^2 - uold(j-1)^2);
end              
subplot(223), hold on, plot(y,u,'o')
title('Lax-Friedrichs scheme','fontsize',14)

u = uinit;  fplus = zeros(1,J); fminus = fplus;	
for n = 1:ntime		% Engquist-Osher scheme (9.82)
  uold = u;
  for j = 1:J
    if uold(j) > 0,	 fplus(j) = 0.5*uold(j)^2; fminus(j) = 0;
    else,		 fplus(j) = 0; fminus(j) = 0.5*uold(j)^2;
    end
  end
  if uleft > 0,	 	fplusleft = 0.5*uleft^2;
  else,		 	fplusleft = 0; 
  end
  if uright > 0, 	fminusright = 0;
  else,		 	fminusright = 0.5*uright^2;
  end
  j = 1; 
  u(j) = uold(j) - lambda*(fplus(j) + fminus(j+1) - fplusleft - fminus(j));  
  for j = 2:J-1   
    u(j) = uold(j) - lambda*(fplus(j) + fminus(j+1) - fplus(j-1) - fminus(j));
  end
  j = J;
  u(j) = uold(j) - lambda*(fplus(j) + fminusright - fplus(j-1) - fminus(j));    
end              
subplot(224), hold on, plot(y,u,'o')
title('Engquist-Osher scheme','fontsize',14)

	% Exact solution
zz = -1:0.01:1; uexact = zeros(size(zz)); t = n*dt;
if geval == 1
  for j = 1:length(zz),
    if zz(j) < t/2,  	uexact(j) = 1;
    else,		uexact(j) = 0;
    end
  end  
else
  for j = 1:length(zz),
    if zz(j) < - t,	uexact(j) = -1;
    elseif zz(j) > t,	uexact(j) = 1;
    else,		uexact(j) = zz(j)/t;
    end
  end  
end
for k = 1:4, subplot(2,2,k), plot(zz,uexact), end
%print solution -depsc;	% Figure is saved in solution.eps for printing

