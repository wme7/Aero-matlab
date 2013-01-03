% Staggered Mach-uniform scheme for one-dimensional Euler equations
% Flux limiting (MUSCL) or first order upwind
% SHK Runge-Kutta time stepping
% Pressure correction method
% Dimensionless formulation; dimensional quantities have postfix dim

% Theory in Section 14.6 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 14.14 -- 14.19 in the book

% Programs called:   Riemann, eulerstep_stag_unif

clear all
global PRL  CRL  MACHLEFT  gamma	% Used in Riemann.m to call f.m

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 49;		% Number of grid cells
limtype = 2;	% Enter 1 for minmod 
		%   or  2 for van Albada 
		%   or  3 for superbee
		%   or  4 for PL limiter (Sect. 4.8) 
		%   or  5 for Chakravarty/Osher limiter (Kr\"oner (1997) p. 109)
		%   or something else for upwind (no MUSCL)
shk = 1;	% Enter 1 for SHK Runge-Kutta method or something else for
		% explicit Euler

problem_specification

	% Dimensional quantities		
pleftdim    = pleft;   prightdim = pright; rholeftdim = rholeft; 
rhorightdim = rhoright; uleftdim = uleft ;  urightdim = uright; 
tenddim = tend; lambdadim = lambda; 

	% Physical units
rhoref = 0.5*(rholeftdim + rhorightdim); pref = 0.5*(pleftdim + prightdim);
uref = sqrt(gamma*pref/rhoref); mrefsq = 1;
tref = 1/uref; lambda = lambdadim/tref;

	% Dimensionless initial conditions
pleft = (pleftdim -pref)/(rhoref*uref^2);
pright = (prightdim -pref)/(rhoref*uref^2);
rholeft = rholeftdim/rhoref; rhoright = rhorightdim/rhoref;
uleft = uleftdim/uref;  uright = urightdim/uref;
tend = tenddim/tref;
 
h = 1/(J-1);  				% Cell size  
dt = lambda*h;				% Time step
n = floor(tend/dt);			% Number of time-steps

% 		Definition of grid numbering 
%     x=0    					x=1
% grid |---o---|---o---|---o---  ....  --|---o---|
%      1   1   2   2   3   3           J-1 J-1  J

xcenter = h*[1:J-1]' - h/2;		% Location of cell centers
xwall = h*[0:J-1]';			% Location of cell boundaries

pold = zeros(size(xcenter));		% Preallocation of pressure, 
rhoold = pold;				%       density

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5 - 5*eps,    pold(j) = pleft;  rhoold(j) = rholeft;
  elseif abs(xcenter(j) - 0.5) < 10*eps
    pold(j) = 0.5*(pleft + pright);
    rhoold(j) = 0.5*(rholeft + rhoright);
  else,    pold(j) = pright;  rhoold(j) = rhoright;
  end
end

uold = zeros(size(xwall));		% Preallocation of velocity, 
mold = uold; um = uold;			%	momentum

for j = 1:length(xwall)			% Initial conditions
  if xwall(j) < 0.5 - 5*eps,    	  uold(j) = uleft;  
  elseif abs(xwall(j) - 0.5) < 10*eps,    uold(j) = 0.5*(uleft + uright);  
  else,    				  uold(j) = uright;
  end
end
for j = 1:length(xwall)
  if xwall(j) < 0.5 - 5*eps,    mold(j) = uleft*rholeft;  
  elseif abs(xwall(j) - 0.5) < 10*eps
  			        mold(j) = 0.5*(uleft*rholeft+uright*rhoright);	
  else,				mold(j) = uright*rhoright;
  end
end

	% Preallocations
pnew = pold; dp = pnew; rhonew = dp; rhoL = dp; pL = dp;
unew = uold; mnew = unew; rhowall = unew;

t = 0;
for i = 1:n
  rhostar = rhoold; mstar = mold; ustar = uold; pstar = pold;
  if shk == 1			% SHK Runge-Kutta time stepping
    rkalpha = 0.25;  eulerstep_stag_unif
    rkalpha = 1/3;   eulerstep_stag_unif
    rkalpha = 0.5;   eulerstep_stag_unif
    rkalpha = 1;     eulerstep_stag_unif
  else				% Explicit Euler time stepping
    rkalpha = 1;     eulerstep_stag_unif
  end
	% Preallocation of J-1*J-1 tridiagonal matrix A
  A = spdiags([ ones(J-1,1)  ones(J-1,1)  ones(J-1,1)], [-1 0 1]',J-1,J-1);
  rhs = zeros(J-1,1); 	% Preallocation of right-hand side 
  
  allamsq = 0.5*allam^2;	% allam is lambda*Runge-Kutta coefficient alpha
  alm = allamsq*mrefsq;
  pL = extrap(pold, limtype);

	% Generation of pressure correction matrix

  A(1,2) = -allamsq/rhowall(2) - (alm/rhowall(2))*(pL(2)+(gamma-1)*pold(1));
  A(1,1) = mrefsq - A(1,2);
  for j = 2:(J-2)
    A(j,j-1) = -allamsq/rhowall(j)- (alm/rhowall(j))*(pL(j)+(gamma-1)*pold(j));
    A(j,j+1) = -allamsq/rhowall(j+1)...
               - (alm/rhowall(j+1))*(pL(j+1)+(gamma-1)*pold(j));
    A(j,j)   = - A(j,j-1) - A(j,j+1) + mrefsq;
  end
  j = J-1;	% Application of pressure boundary conditions
  A(j,j-1) = - allamsq/rhowall(j)- (alm/rhowall(j))*(pL(j)+(gamma-1)*pold(j));
  A(j,j)   = mrefsq - A(j,j-1) + 2*allamsq/rhowall(j+1)...
             + 2*(alm/rhowall(j+1))*(pright+(gamma-1)*pold(j));
	       
	% Right-hand side for pressure correction equation
	
  rhs(1) = - mrefsq*allam*(ustar(2)*pL(1) - uleft*pleft + (gamma-1)*...
             pstar(1)*(ustar(2) - uleft));
  rhs(1) = rhs(1)  - allam*(ustar(2) - uleft);
  for j = 2:(J-2)
    rhs(j) =  - mrefsq*allam*(ustar(j+1)*pL(j) - ustar(j)*pL(j-1) + ...
  	     (gamma-1)*pstar(j)*(ustar(j+1) - ustar(j)));
    rhs(j) = rhs(j) - allam*(ustar(j+1) - ustar(j));
  end
  j = J-1;
  rhs(j) = - mrefsq*allam*(ustar(j+1)*pright - ustar(j)*pL(j-1) + ...
           (gamma-1)*pstar(j)*(ustar(j+1) - ustar(j)));
  rhs(j) = rhs(j) - allam*(ustar(j+1) - ustar(j));

  dp = A\rhs;  pstar = pold + dp;

  for j = 2:(J-1)
    mstar(j) = mstar(j) - 0.5*allam*(dp(j) - dp(j-1));
  end
  mstar(J) = mstar(J) + allam*dp(J-1);
  ustar = mstar./rhowall;
  pnew = pstar; rhonew = rhostar; mnew = mstar; unew = ustar;
  pold = pnew; rhoold = rhonew; uold = unew; mold = mnew; 
end

	% Dimensional quantities for output
pnewdim = pref*ones(size(pnew)) + rhoref*uref^2*pnew;
rhonewdim = rhoref*rhonew; unewdim = uref*unew; tdim = n*dt*tref;
mnewdim = mnew*uref*rhoref;

entropy = log(pnewdim./(rhonewdim.^gamma));
soundspeed = sqrt(gamma.*pnewdim./rhonewdim);
cc = zeros(size(xwall));     % Preallocation of sound speed at cell boundaries 
cc(1) = sqrt(gamma.*pleftdim./rholeftdim);
cc(length(xwall)) = sqrt(gamma.*prightdim./rhorightdim);
for j = 2:length(xwall) - 1
  cc(j) = 0.5*(soundspeed(j-1) + soundspeed(j));
end
mach  = unewdim./cc;

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14)
	plot(xcenter,rhonewdim,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14)
	plot(xwall,unewdim,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14)
	plot(xcenter,pnewdim,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xwall,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6),hold on, axis('off'), title('Staggered scheme')
text(0,0.9,'Pressure correction scheme 2')
if shk == 1,   text(0,0.75,'SHK Runge-Kutta')
else, 	       text(0,0.75,'explicit Euler')
end
if limtype == 0,   s1 = 'No MUSCL';
elseif limtype == 1,  s1 = 'minmod limiter';
elseif limtype == 2,  s1 = 'van Albada limiter';
elseif limtype == 3,  s1 = 'superbee limiter';
elseif limtype == 4,  s1 = 'PL limiter';
elseif limtype == 5,  s1 = 'Chakravarty-Osher limiter';
else,  s1 = 'limtype wrong';
end
text(0,0.6,s1)
text(0,0.45,['t = ', num2str(tdim,3)])
text(0,0.3,['lambda = ', num2str(lambdadim),' J = ',num2str(J)])
%klok = clock;
%s = [date,' ',num2str(klok(4)),':',num2str(klok(5))];
%text(0,0,s)

% Change to dimensional quantities before call of Riemann
pright = prightdim; pleft = pleftdim;
uright = urightdim; uleft = uleftdim;
rhoright = rhorightdim; rholeft =  rholeftdim;
tend = tdim;
Riemann

