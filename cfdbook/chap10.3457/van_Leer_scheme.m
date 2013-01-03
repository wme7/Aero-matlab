% Van Leer scheme for one-dimensional Euler equations

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.5 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.15, 10.16 in the book

% Functions called: f, problem_specification, Riemann 

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
totenew = press;			%	total energy

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

% Preallocation of split fluxes
plus1 = c; minus1 = c; plus2 = c; minus2 = c; plus3 = c; minus3 = c;
% Preallocation of van Leer fluxes
flux1 = zeros(J-1,1); flux2 = flux1; flux3 = flux1;

t = 0;
for i = 1:n,  t = t + dt;
  Eflux1 = rhoold.*uold;			% Eflux1,2,3 is Euler flux
  Eflux2 = Eflux1.*uold + (1/gamma)*rhoold.*c.^2;
  Eflux3 = 0.5*Eflux1.*uold.^2 + gammab*Eflux1.*c.^2;  
  
    for j = 1:J
    if mach(j) > 1
      plus1(j) = Eflux1(j);  plus2(j) = Eflux2(j);  plus3(j) = Eflux3(j); 
    elseif mach(j) < -1
      plus1(j) = 0;         plus2(j) = 0;         plus3(j) = 0; 
    else
      plus1(j) = 0.25*rhoold(j)*c(j)*(1 + mach(j))^2;
      plus2(j) = plus1(j)*c(j)*(2 + gam1*mach(j))/gamma;
      plus3(j) = (plus2(j)^2/plus1(j))*gamma^2*0.5*gammab/(gamma+1);
    end
  end
  minus1 = Eflux1 - plus1; minus2 = Eflux2 - plus2; minus3 = Eflux3 - plus3;    

  for j = 1:J-1					% van Leer fluxes
    flux1(j) = plus1(j) + minus1(j+1);
    flux2(j) = plus2(j) + minus2(j+1);
    flux3(j) = plus3(j) + minus3(j+1);
  end
     
	% Update of state variables
  rhonew(1) = rholeft;		  rhonew(J)  = rhoright; 
  mnew(1) = rholeft*uleft; 	  mnew(J)    = rhoright*uright;
  totenew(1) = totenleft; 	  totenew(J) = totenright;
  for j = 2:J-1
    rhonew(j)  = rhoold(j)   - lambda*(flux1(j) - flux1(j-1));
    mnew(j)    = mold(j)     - lambda*(flux2(j) - flux2(j-1));
    totenew(j) = totenold(j) - lambda*(flux3(j) - flux3(j-1));
  end

  if bouncon ~= 0
    if bouncon == 1
      rhowall = rhoold(J); uwall = 0;
      pwall = (gamma-1)*(totenold(J) - 0.5*mold(J)^2/rhoold(J)); % pwall = p(J)  
      totenwall = pwall/(gamma-1) + 0.5*rhowall*uwall^2;
    else
      rhowall = rhoold(J); uwall = - mold(J)/rhoold(J);
      pwall = (gamma-1)*(totenold(J) - 0.5*mold(J)^2/rhoold(J)); % pwall = p(J)         
      totenwall = pwall/(gamma-1) + 0.5*rhowall*uwall^2;
    end
    wallflux1 = rhowall*uwall;
    pwall = gam1*(totenwall - 0.5*wallflux1*uwall);
    wallflux2 = wallflux1.*uwall + pwall;
    wallflux3 = 0.5*wallflux1*uwall^2 + gammab*gamma*uwall*pwall;
    cwall = sqrt(gamma*pwall/rhowall); machwall = uwall/cwall;  
    if machwall > 1
      plus1wall = wallflux1; plus2wall = wallflux2; plus3wall = wallflux3 ; 
    elseif machwall < -1
      plus1wall = 0;      plus2wall = 0;      plus3wall = 0; 
    else
      plus1wall = 0.25*rhowall*cwall*(1 + machwall)^2;
      plus2wall = plus1wall*cwall*(2 + gam1*machwall)/gamma;
      plus3wall = (plus2wall^2/plus1wall)*gamma^2*0.5*gammab/(gamma+1);
    end
    minus1wall = wallflux1 - plus1wall; 
    minus2wall = wallflux2 - plus2wall; minus3wall = wallflux3 - plus3wall;
        
    flux1wall = plus1(J) + minus1wall;			% van Leer fluxes
    flux2wall = plus2(J) + minus2wall;
    flux3wall = plus3(J) + minus3wall;
    
    rhonew(J)  = rhoold(J)   - lambda*(flux1wall - flux1(J-1));
    mnew(J)    = mold(J)     - lambda*(flux2wall - flux2(J-1));
    totenew(J) = totenold(J) - lambda*(flux3wall - flux3(J-1));
  end
 	% Update of state variables 
  uold = mnew./rhonew; press = gam1*(totenew - 0.5*mnew.*uold);
  rhoold = rhonew; totenold = totenew; mold = mnew;
  
  c = sqrt(gamma*press./rhoold); mach = uold./c;
end

entropy = log(press./rhoold.^gamma);

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6),axis('off'),hold on,title('Van Leer scheme','fontsize',14)

Riemann		% Plot exact solution

