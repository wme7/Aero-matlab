% Roe scheme with MUSCL for one-dimensional Euler equations
% Runge-Kutta time stepping

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.28 -- 10.31 in the book

% Functions called: extrap, roeeulerstep, problem_specification, Riemann 

clear all
global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  epsi  invariants... 
	lambda		% lambda = dt/dx

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 48;		% Number of grid cells
invariants = 0;	% Enter 0 for extrapolation of primitive variables or
		%   something else for extrapolation of Riemann invariants 
limtype = 2;	% Enter 1 for minmod or 2 for van Albada or something else
		%	if MUSCL not wanted
epsi = 0.0;	% Parameter in Harten's sonic entropy fix, see Sect. 10.3
		% ....................End of input........................

gammab = 1/(gamma - 1); gam1 = gamma-1; 
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

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5, press(j) = pleft;  rhoold(j) = rholeft; uold(j) = uleft;  
  else,		      press(j) = pright; rhoold(j) = rhoright; uold(j) = uright;
  end
end

	% Initialization of cell center variables
eold = gammab*press./rhoold;		% e = internal energy
enrhoold = gammab*press;		% enrho = rho*e
mold = uold.*rhoold;			% m = momentum
totenold = enrhoold + 0.5*mold.*uold;	% toten = rho*(total energy)
totenleft = totenold(1); totenright = totenold(J);
htot = gamma*eold + 0.5*uold.*uold;	% htot = total enthalpy

	% Preallocation of cell center variables
rhonew = press; mach = press; 
 mnew = press; enrhonew = press; totenew = press; cnew = press;
Z1 = press; Z2 = Z1; Z3 = Z1;		% Z1, Z2, Z3: Riemann variables
	% Preallocation of cell edge variables
roeflux1 = zeros(J-1,1); roeflux2 = roeflux1; roeflux3 = roeflux1; % Roe fluxes
hav = roeflux1; uav = hav; cav = hav; dd = hav;		% Roe averages
f1av = hav; f2av = hav; f3av = hav;			% Flux averages
delrho = hav; delm = hav; deltoten = hav;	% State vector differences
alambda1 = hav; alambda2 = hav; alambda3 = hav;	% Eigenvalues and coefficients
alpha1 = hav; alpha2 = hav; alpha3 = hav;
Z1L = zeros(J-1,1); Z2L = Z1L; Z3L = Z1L;		% Extrapolated states
Z1R = Z1L; Z2R = Z1L; Z3R = Z1L;
U1L = Z1L; U1R = Z1L; U2L = Z1L; U2R= Z1L; U3R = Z1L; U3L = Z1L;
 
t = 0;
for i = 1:n,   t = t + dt;
  rhostar = rhoold; mstar = mold; totenstar = totenold;
  rkalpha = 0.25; roeeulerstep
  rkalpha = 1/3;  roeeulerstep
  rkalpha = 0.5;  roeeulerstep
  rkalpha = 1;    roeeulerstep
  rhonew = rhostar; 	mnew = mstar; 	totenew = totenstar;  
  uold = mnew./rhonew;  press = gam1*(totenew - 0.5*mnew.*uold);
  eold = gammab*press./rhonew;	  	htot = (totenew + press)./rhonew;
  rhoold = rhonew;	mold = mnew;	totenold = totenew;
  cnew = sqrt(gamma*gam1*eold); 	mach = uold./cnew;
end

entropy = log(press./rhoold.^gamma);
figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o')
subplot(2,3,6), axis('off'), hold on, title('Roe scheme','fontsize',14)
text(0,0.9,['lambda = ', num2str(lambda),'  t = ',num2str(n*dt)])
if invariants == 0,  s1 = 'Primitive extrapolation';
else,  s1 = 'Riemann extrapolation';
end
text(0,0.75,s1),	text(0,0.6,'SHK Runge-Kutta')
if limtype == 0,  	s1 = 'No MUSCL';
elseif limtype == 1,  	s1 = 'MUSCL, minmod limiter';
else,  			s1 = 'MUSCL, van Albada limiter';
end
text(0,0.45,s1),	text(0,0.3,['epsilon = ',num2str(epsi)])

Riemann		% Plot exact solution
