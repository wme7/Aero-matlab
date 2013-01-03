% AUSM (Liou-Steffen) scheme for one-dimensional Euler equations

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.5 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.17, 10.20 in the book

% Function called: f
% Programs called: problem_specification, Riemann

clear all

global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  lambda		% lambda = dt/dx

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 48;		% Number of grid cells
bouncon = 0;	% bouncon chooses outflow boundary conditions
		% = 0: Nothing happens: infinite domain
		% = 1: Solid wall at x = 1 with direct prescription of uwall
		% = 2: Solid wall at x = 1 with reflection b.c.
		% ....................End of input........................

gammab = 1/(gamma - 1); gam1 = gamma-1; gamgam = gamma^gamma;
problem_specification	
		
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
totenew = press; enthalpy = press;	%	total energy and enthalpy

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5, press(j) = pleft; rhoold(j) = rholeft;  uold(j) = uleft;  
  else,  	     press(j) = pright;  rhoold(j) = rhoright; uold(j) = uright;
  end
end

	% Initialization of cell center variables
totenold = rhoold.*(0.5*uold.*uold + gammab*press./rhoold); % Total energy rho*E
totenleft = totenold(1); totenright = totenold(J);
mold = rhoold.*uold;					    % Momentum m
c = sqrt(gamma*press./rhoold);				    % Sound speed 
mach = uold./c;						    % Mach number
enthalpy =  0.5*uold.*uold + gammab*c.^2;		    % Enthalpy

machplus = mach; machminus = mach;		% Preallocation of  
presplus = mach; presminus = mach;		%     	split fluxes 
flux1 = zeros(J-1,1); flux2 = flux1;		% 	and Liou-Steffen 
flux3 = flux1; machhalf = flux1;		%	fluxes
m1 = flux1; m2 = flux1;

t = 0;
for i = 1:n,  t = t + dt;    
  for j = 1:J
    if mach(j) > 1
      machplus(j) = mach(j);  machminus(j) = 0;   
      presplus(j) = press(j); presminus(j) = 0;
    elseif mach(j) < -1
      machplus(j) = 0; machminus(j) = mach(j);   
      presplus(j) = 0; presminus(j) = press(j);
    else
      machplus(j)  =  0.25*(mach(j) + 1)^2;
      machminus(j) = -0.25*(mach(j) - 1)^2;  
      presplus(j)  =  0.5*press(j)*(1 + mach(j));
      presminus(j) =  0.5*press(j)*(1 - mach(j));
    end
  end

	% Liou-Steffen fluxes
  for j = 1:J-1,  machhalf(j) = machplus(j) + machminus(j+1); end
  m1 = 0.5*(machhalf + abs(machhalf));  m2 = 0.5*(machhalf - abs(machhalf));
  rhoc = rhoold.*c;
  for j = 1:J-1
    flux1(j) = m1(j)*rhoc(j) + m2(j)*rhoc(j+1);    
    flux2(j) = m1(j)*rhoc(j)*uold(j) +  m2(j)*rhoc(j+1)*uold(j+1) +...
      presplus(j) + presminus(j+1); 
    flux3(j) =  m1(j)*rhoc(j)*enthalpy(j)  + m2(j)*rhoc(j+1)*enthalpy(j+1);   
  end
     
	% Update of state variables
  rhonew(1)  = rholeft;  	rhonew(J)  = rhoright; 
  mnew(1)    = rholeft*uleft;	mnew(J)    = rhoright*uright;
  totenew(1) = totenleft; 	totenew(J) = totenright;
  for j = 2:J-1
    rhonew(j)  = rhoold(j)   - lambda*(flux1(j) - flux1(j-1));
    mnew(j)    = mold(j)     - lambda*(flux2(j) - flux2(j-1));
    totenew(j) = totenold(j) - lambda*(flux3(j) - flux3(j-1));
  end
  uold   = mnew./rhonew;   press = gam1*(totenew - 0.5*mnew.*uold);
  rhoold = rhonew; 	totenold = totenew; 	mold = mnew;
  c = sqrt(gamma*press./rhoold); 	mach = uold./c;
  enthalpy =  0.5*uold.*uold + gammab*c.^2;
end
entropy = log(press./rhoold.^gamma);

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6),axis('off'),title('Liou-Steffen scheme','fontsize',14)

Riemann		% Plot exact solution
