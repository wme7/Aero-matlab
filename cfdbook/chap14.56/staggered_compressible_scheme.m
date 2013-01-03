% Staggered compressible scheme for one-dimensional Euler equations
% Flux limiting (MUSCL) or first order upwind
% SHK Runge-Kutta time stepping

% Theory in Section 14.5 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 14.2 -- 14.7 in the book

% Programs called:   problem_specification, Riemann, eulerstep_stag_com

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
theta = 0;	% Governs upwind bias in mass conservation equation.
		% Enter 0 for full upwind or 1 for central
shk = 1;	% Enter 1 for SHK Runge-Kutta method or something else for
		% explicit Euler
		%...................End of input.......................
		
problem_specification	

h = 1/(J-1);  				% Cell size  
dt = lambda*h;				% Time step
n = floor(tend/dt);			% Number of time-steps

% 		Definition of grid numbering 
%     x=0					x=1    
% grid |---o---|---o---|---o---  ....  --|---o---|
%      1   1   2   2   3   3           J-1 J-1  J

xcenter = h*[1:J-1] - h/2;		% Location of cell centers
xwall = h*[0:J-1];			% Location of cell boundaries

toll = 10*eps;	% Used to specify precisely the initial conditions for   
shift = 0;	% quantities in a gridpoint that is located at the jump in the 
		% initial conditions for the Riemann problems

pold = zeros(size(xcenter));		% Preallocation of pressure, 
rhoold = zeros(size(xcenter));		%       density

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5 + shift - toll
    pold(j) = pleft;  rhoold(j) = rholeft;
  elseif xcenter(j) > 0.5 + shift + toll
    pold(j) = pright;  rhoold(j) = rhoright;
  else
    pold(j) = 0.5*(pright+pleft);  rhoold(j) =  0.5*(rhoright+rholeft);
  end
end

uold = zeros(size(xwall));		% Preallocation of velocity, 
mold = zeros(size(xwall));		%	momentum

for j = 1:length(xwall)			% Initial conditions
  if xwall(j) < 0.5  + shift- toll
    uold(j) = uleft;    mold(j) = uleft*rholeft; 
  elseif xwall(j) > 0.5  + shift+ toll
    uold(j) = uright;    mold(j) = uright*rhoright;
  else
    uold(j) = 0.5*(uright+uleft);
    mold(j) = 0.25*(uright+uleft)*(rhoright+rholeft);
%    uold(j) = 2*mold(j)/(rhoright+rholeft);
  end
end

etotold = zeros(size(xcenter));		% Preallocation of total energy
eint = (1/(gamma-1))*pold;	% eint = rho * internal energy

for j = 1:length(xcenter)		% Initial conditions
    etotold(j) = eint(j) + 0.25*(mold(j)*uold(j) + mold(j+1)*uold(j+1));
end

	% Preallocations
pnew = zeros(size(xcenter));     rhonew = zeros(size(xcenter));
etotnew = zeros(size(xcenter)); 
unew = zeros(size(xwall));         mnew = zeros(size(xwall));
totenthL = zeros(size(xcenter)); totenthalpy = zeros(size(xcenter));
umL = zeros(size(xwall)); um = zeros(size(xwall)); rhoL = zeros(size(xcenter));
rhostar = zeros(size(xcenter));                   mstar = zeros(size(xwall));
etotstar = zeros(size(xcenter));                  ustar = zeros(size(xwall));

t = 0;
for i = 1:n,  t = t + dt;
  rhostar = rhoold; mstar = mold; etotstar = etotold; ustar = uold;
  pstar = pold;
  if shk == 1			% SHK Runge-Kutta time stepping
    rkalpha = 0.25;  eulerstep_stag_com
    rkalpha = 1/3;   eulerstep_stag_com
    rkalpha = 0.5;   eulerstep_stag_com
    rkalpha = 1;     eulerstep_stag_com
  else				% Explicit Euler time stepping
    rkalpha = 1;     eulerstep_stag_com
  end
  rhonew = rhostar; mnew = mstar; etotnew = etotstar; unew = ustar;
  pnew   = pstar;   pold = pnew;   rhoold = rhonew;   uold = unew; 
  mold   = mnew; etotold = etotnew;
end 

entropy = log(pnew./(rhonew.^gamma)); soundspeed = sqrt(gamma.*pnew./rhonew);
cc = zeros(size(xwall));		% Sound speed at cell boundaries
cc(1)             = sqrt(gamma.*pleft./rholeft);
cc(length(xwall)) = sqrt(gamma.*pright./rhoright);
for j = 2:length(xwall) - 1
  cc(j) = 0.5*(soundspeed(j-1) + soundspeed(j));
end
mach  = unew./cc;

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xwall,unew,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,pnew,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xwall,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6), axis('off'), hold on
	title('Staggered scheme','fontsize',14)
text(0,0.9,'Conservation form')
if shk ==1,  text(0,0.75,'SHK Runge-Kutta')
else,  text(0,0.75,'explicit Euler')
end
if limtype == 0,  	s1 = 'No MUSCL';
elseif limtype == 1,  	s1 = 'Minmod limiter';
elseif limtype == 2,  	s1 = 'van Albada limiter';
elseif limtype == 3,  	s1 = 'Superbee limiter';
else,			s1 = 'wrong limtype'
end
text(0,0.6,s1)
text(0,0.45,[' t = ', num2str(n*dt,3)])
text(0,0.3,['lambda = ', num2str(lambda),'  J = ', num2str(J)])
%klok = clock;
%s = [date,' ',num2str(klok(4)),':',num2str(klok(5))];
%text(0,0,s)

Riemann		% Plot exact solution

