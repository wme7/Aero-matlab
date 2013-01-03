% Jameson-Schmidt-Turkel scheme for one-dimensional Euler equations

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.7 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.22 -- 10.26 in the book

% Functions called: f, problem_specification, Riemann 

global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  lambda		% lambda = dt/dx

			% ...........Input................................
gamma = 1.4; 		% Ratio of specific heats
J = 48;			% Number of grid cells
k2 = 8; 		% Jameson's artificial diffusion coefficients
k4 = 1/32; 		% according to equation (10.98), for example.
alpjam = 2;		% alpha in equation (10.98)
			% ..........End of input..........................
			
gammab = 1/(gamma - 1); gam1 = gamma-1; 
problem_specification
	
alplam1 = lambda/3; 	% Runge-Kutta time stepping 
alplam2 = 4*lambda/15;	%   coefficients, equation (10.101)
alplam3 = 5*lambda/9;	%     alpha*lambda
		
h = 1/J;  				% Cell size
dt = lambda*h;				% Time step
n = floor(tend/dt);			% Number of time-steps

% 		Definition of grid numbering 
%       x=0    					 x=1
% grid   |---o---|---o---|---o---  ...  --|---o---|
%            1   1   2   2   3           J-1  J  

xcenter = h*[1:J] - h/2;		% Location of cell centers

press = zeros(size(xcenter));		% Preallocation of pressure,  
rhoold = press; uold = press;		%       density and velocity
rhonew = press; mnew = press;		%       momentum 
totenew = press;			%	total energy

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5, press(j) = pleft; rhoold(j) = rholeft;  uold(j) = uleft;  
  else,  	     press(j) = pright;  rhoold(j) = rhoright; uold(j) = uright;
  end
end

totenold = rhoold.*(0.5*uold.*uold + gammab*press./rhoold); % Total energy rho*E
mold = rhoold.*uold;					    % Momentum m

	% Preallocations
rjac = zeros(J-1,1)';			% Spectral radius of Jacobian 
eps2 = rjac; eps4 = rjac;	% Jameson's coefficients above equation (10.98)
d1 = rjac; d2 = rjac;d3 = rjac;		% Diffusive fluxes at cell faces
	% Central and diffusive fluxes in cell centers
cenflux1 = press; cenflux2 = cenflux1; cenflux3 = cenflux1;
difflux1 = press; difflux2 = difflux1; difflux3 = difflux1;
	% Shock sensor in cell centers and at cell faces
nu = press; nubar = rjac;

c = sqrt(gamma*press./rhoold);		% Sound speed 

t = 0;
for i = 1:n,  t = t + dt;
  
  Eflux1 = rhoold.*uold;		% Eflux1,2,3 are Euler fluxes
  Eflux2 = Eflux1.*uold + press;
  rhoent = rhoold.*(gammab*c.^2 + 0.5*uold.^2);	% rhoent = rho*enthalpy
  Eflux3 = uold.*rhoent;
  
  cenflux1(1) = 0; cenflux2(1) = 0; cenflux3(1) = 0;
  for j = 2:J-1
    cenflux1(j) = 0.5*(Eflux1(j+1) - Eflux1(j-1));
    cenflux2(j) = 0.5*(Eflux2(j+1) - Eflux2(j-1));
    cenflux3(j) = 0.5*(Eflux3(j+1) - Eflux3(j-1));
  end
  cenflux1(J) = 0; cenflux2(J) = 0; cenflux3(J) = 0;
  
  for j = 1:J-1			% Spectral radius of Jacobian at cell faces
    rjac(j) = 0.5*(abs(uold(j)) + abs(uold(j+1)) + c(j) + c(j+1));
  end
  
	% Shock sensor in cell centers
  nu(1) = abs((press(2) - press(1))/(press(2)+3*press(1)));
  for j = 2:J-1
    nu(j) = press(j+1)-2*press(j)+press(j-1);
    nu(j) = nu(j)/(press(j+1)+2*press(j)+press(j-1)); nu(j) = abs(nu(j));
  end
  nu(J) = abs((press(J-1)-press(J))/(press(J-1)+3*press(J)));
  
  nubar(1) = max(nu(3),max(nu(2),nu(1)));	% Shock sensor at cell faces
  for j = 2:J-2
    nubar(j) = max(nu(j-1),max(nu(j),max(nu(j+1),nu(j+2))));
  end
  nubar(J-1) = max(nu(J-2),max(nu(J-1), nu(J)));
  
  for j = 1:J-1			% Artificial viscosities
    eps2(j) = min(0.5, k2*nubar(j)); eps4(j) = max(0, k4 - alpjam*nubar(j));
  end
  
	% Diffusive fluxes at cell faces
  d1(1) = rjac(1)*(eps2(1)*(rhoold(2) - rhoold(1))...
     - eps4(1)*(rhoold(3) -3*rhoold(2)+ 2*rhoold(1)));
  d2(1) = rjac(1)*(eps2(1)*(mold(2) - mold(1))...
     - eps4(1)*(mold(3) -3*mold(2)+ 2*mold(1)));
  d3(1) = rjac(1)*(eps2(1)*(rhoent(2) - rhoent(1))...
     - eps4(1)*(rhoent(3) -3*rhoent(2)+ 2*rhoent(1)));
  for j = 2:J-2
    d1(j) = rjac(j)*(eps2(j)*(rhoold(j+1) - rhoold(j))...
      - eps4(1)*(rhoold(j+2) -3*rhoold(j+1)+ 3*rhoold(j) -rhoold(j-1)));
    d2(j) = rjac(j)*(eps2(j)*(mold(j+1) - mold(j))...
     - eps4(j)*(mold(j+2) -3*mold(j+1)+ 3*mold(j) - mold(j-1)));
    d3(j) = rjac(j)*(eps2(j)*(rhoent(j+1) - rhoent(j))...
     - eps4(j)*(rhoent(j+2) -3*rhoent(j+1)+ 3*rhoent(j) - rhoent(j-1)));
  end
  d1(J-1) = rjac(J-1)*(eps2(J-1)*(rhoold(J) - rhoold(J-1))...
     - eps4(J-1)*( -2*rhoold(J)+ 3*rhoold(J-1) - rhoold(J-2)));
  d2(J-1) = rjac(J-1)*(eps2(J-1)*(mold(J) - mold(J-1))...
     - eps4(J-1)*(-2*mold(J)+ 3*mold(J-1)-mold(j-2)));     
  d3(J-1) = rjac(J-1)*(eps2(J-1)*(rhoent(J) - rhoent(J-1))...
     - eps4(J-1)*( -2*rhoent(J)+ 3*rhoent(J-1)- rhoent(J-2)));
 
	% Diffusive fluxes at cell centers
    difflux1(1) = 0; difflux2(1) = 0; difflux3(1) = 0;
  for j = 2:J-1, difflux1(j) = d1(j) - d1(j-1);
                 difflux2(j) = d2(j) - d2(j-1); difflux3(j) = d3(j) - d3(j-1); 
  end
    difflux1(J) = 0; difflux2(J) = 0; difflux3(J) = 0;
 
	% First Runge-Kutta stage
  rho1 = rhoold - alplam1*(cenflux1 -difflux1);
  Eflux1 = mold - alplam1*(cenflux2 -difflux2);
  toten1 = totenold - alplam1*(cenflux3 -difflux3);  
	% Eflux1,2,3 are Euler fluxes
  u1 = Eflux1./rho1; Eflux2 = Eflux1.*u1 + gam1*(toten1 - 0.5*rho1.*u1.^2);
  Eflux3 = u1.*(gamma*toten1 - 0.5*gam1*rho1.*u1.^2);
  
  cenflux1(1) = 0; cenflux2(1) = 0; cenflux3(1) = 0;
  for j = 2:J-1
    cenflux1(j) = 0.5*(Eflux1(j+1) - Eflux1(j-1));
    cenflux2(j) = 0.5*(Eflux2(j+1) - Eflux2(j-1));
    cenflux3(j) = 0.5*(Eflux3(j+1) - Eflux3(j-1));
  end
  cenflux1(J) = 0; cenflux2(J) = 0; cenflux3(J) = 0;
  
	% Second Runge-Kutta stage
  rho1 = rhoold - alplam2*(cenflux1 -difflux1);
  Eflux1 = mold - alplam2*(cenflux2 -difflux2);
  toten1 = totenold - alplam2*(cenflux3 -difflux3);
	% Eflux1,2,3 are Euler fluxes
  u1 = Eflux1./rho1; Eflux2 = Eflux1.*u1 + gam1*(toten1 - 0.5*rho1.*u1.^2);
  Eflux3 = u1.*(gamma*toten1 - 0.5*gam1*rho1.*u1.^2);
  
  cenflux1(1) = 0; cenflux2(1) = 0; cenflux3(1) = 0;
  for j = 2:J-1
    cenflux1(j) = 0.5*(Eflux1(j+1) - Eflux1(j-1));
    cenflux2(j) = 0.5*(Eflux2(j+1) - Eflux2(j-1));
    cenflux3(j) = 0.5*(Eflux3(j+1) - Eflux3(j-1));
  end
  cenflux1(J) = 0; cenflux2(J) = 0; cenflux3(J) = 0;
  
	% Third Runge-Kutta stage
  rho1 = rhoold - alplam3*(cenflux1 -difflux1);
  Eflux1 = mold - alplam3*(cenflux2 -difflux2);
  toten1 = totenold - alplam3*(cenflux3 -difflux3);
	% Eflux1,2,3 are Euler fluxes
  u1 = Eflux1./rho1; Eflux2 = Eflux1.*u1 + gam1*(toten1 - 0.5*rho1.*u1.^2);
  Eflux3 = u1.*(gamma*toten1 - 0.5*gam1*rho1.*u1.^2);
  
  cenflux1(1) = 0; cenflux2(1) = 0; cenflux3(1) = 0;
  for j = 2:J-1
    cenflux1(j) = 0.5*(Eflux1(j+1) - Eflux1(j-1));
    cenflux2(j) = 0.5*(Eflux2(j+1) - Eflux2(j-1));
    cenflux3(j) = 0.5*(Eflux3(j+1) - Eflux3(j-1));
  end
  cenflux1(J) = 0; cenflux2(J) = 0; cenflux3(J) = 0;
  
	% Fourth Runge-Kutta stage
  rhonew = rhoold - lambda*(cenflux1 -difflux1);
  mnew = mold - lambda*(cenflux2 -difflux2);
  totenew = totenold - lambda*(cenflux3 -difflux3);
   
  uold = mnew./rhonew; press = gam1*(totenew - 0.5*mnew.*uold);
  rhoold = rhonew; totenold = totenew; mold = mnew;
  c = sqrt(gamma*press./rhoold);
end

mach = uold./c; entropy = log(press./rhoold.^gamma);

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6),axis('off'), hold on, title('JST scheme','fontsize',14)
  text(0,0.9,['lambda = ', num2str(lambda)],'fontsize',14)
  text(0,0.75,['k2 = ',num2str(k2)],'fontsize',14)
  text(0,0.6,['k4 = ',num2str(k4)],'fontsize',14)

Riemann		% Plot exact solution

