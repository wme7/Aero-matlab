% Roe scheme for one-dimensional Euler equations

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.3 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.5 -- 10.10 in the book

% Functions called: f, problem_specification, Riemann 

global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  lambda		% lambda = dt/dx

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 48;		% Number of grid cells
		% ....................End of input........................

gammab = 1/(gamma - 1); gam1 = gamma-1;
problem_specification	
		
h = 1/J;  		% Cell size
dt = lambda*h;		% Time step
n = floor(tend/dt);	% Number of time-steps

% 		Definition of grid numbering 
%       x=0    					 x=1
% grid   |---o---|---o---|---o---  ...  --|---o---|
%            1   1   2   2   3           J-1  J  

xcenter = h*[1:J] - h/2;		% Location of cell centers

press = zeros(size(xcenter));		% Preallocation of pressure,  
rhoold = press; uold = press;		%      density and velocity

for j = 1:length(xcenter)		% Initial conditions
  if xcenter(j) < 0.5, press(j) = pleft; rhoold(j) = rholeft;  uold(j) = uleft; 
  else,		      press(j) = pright; rhoold(j) = rhoright; uold(j) = uright;
  end 
end

	% Initialization of cell center variables
eold = gammab*press./rhoold;		% Internal energy e
enrhoold = gammab*press;		% Internal energy times density
mold = uold.*rhoold;			% Momentum m
totenold = enrhoold + 0.5*mold.*uold;	% Total enery rho*E
totenleft = totenold(1); totenright = totenold(J);
htot = gamma*eold + 0.5*uold.*uold;	% Total enthalpy

rhonew = press; 	% Preallocation of cell center variables
mach = press; mnew = press; enrhonew = press; totenew = press; cnew = press;

			% Preallocation of cell edge variables
roeflux1 = zeros(J-1,1); roeflux2 = roeflux1; roeflux3 = roeflux1; % Roe fluxes
hav = roeflux1; uav = hav; cav = hav; dd = hav;			% Roe averages
f1av = hav; f2av = hav; f3av = hav;				% Flux averages
delrho = hav; delm = hav; deltoten = hav;	% State vector differences
alambda1 = hav; alambda2 = hav; alambda3 = hav;	% Eigenvalues and coefficients
alpha1 = hav; alpha2 = hav; alpha3 = hav;
 
t = 0;
for i = 1:n,  t = t + dt; 
  for j = 1:J-1
    dd(j) = sqrt(rhoold(j+1)/rhoold(j));		% Roe averages
    hav(j) = (htot(j) + dd(j)*htot(j+1))/(1+dd(j));
    uav(j) = (uold(j) + dd(j)*uold(j+1))/(1+dd(j));

    delrho(j) = rhoold(j+1) - rhoold(j);
    delm(j) = mold(j+1) - mold(j);
    deltoten(j) = totenold(j+1) - totenold(j);

    f1av(j) = 0.5*(mold(j) + mold(j+1));
    f2av(j) = 0.5*(mold(j)*mold(j)/rhoold(j) + press(j) +...
      mold(j+1)*mold(j+1)/rhoold(j+1) + press(j+1));
    f3av(j) = 0.5*(mold(j)*htot(j) + mold(j+1)*htot(j+1));
  end
  cav = sqrt(gam1*(hav - 0.5*uav.*uav));  mav = uav./cav;
  alpha1 = 0.25*mav.*(2+gam1*mav).*delrho - 0.5*(1+gam1*mav).*delm./cav +...
           0.5*gam1*deltoten./(cav.^2);
  alpha2 = (1-0.5*gam1*(mav.^2)).*delrho  + gam1*(mav./cav).*delm -...
           gam1*deltoten./(cav.^2);
  alpha3 = -0.25*mav.*(2-gam1*mav).*delrho + 0.5*(1-gam1*mav).*delm./cav +...
           0.5*gam1*deltoten./(cav.^2);
  alambda1 = abs(uav - cav); alambda2 = abs(uav); alambda3 = abs(uav + cav);

  epsi = 0.5;		% Harten's sonic entropy fix
%  epsi = 0;		% No entropy fix
  for j = 1:J-1
    if alambda1(j) < epsi, alambda1(j) = 0.5*(epsi + alambda1(j)^2/epsi); end
    if alambda3(j) < epsi, alambda3(j) = 0.5*(epsi + alambda3(j)^2/epsi); end
  end
  roeflux1 = f1av - 0.5*alambda1.*alpha1...
     - 0.5*alambda2.*alpha2 - 0.5*alambda3.*alpha3;
  roeflux2 = f2av - 0.5*alambda1.*alpha1.*(uav-cav) - ...
     0.5*alambda2.*alpha2.*uav - 0.5*alambda3.*alpha3.*(uav+cav);
  roeflux3 = f3av -...
     0.5*alambda1.*alpha1.*(hav - uav.*cav) -...
     0.25*alambda2.*alpha2.*uav.*uav - 0.5*alambda3.*alpha3.*(hav + uav.*cav);
  rhonew(1) = rholeft; rhonew(J) = rhoright;
  mnew(1) = rholeft*uleft; mnew(J) = rhoright*uright;
  totenew(1) = totenleft; totenew(J) = totenright;
  for j = 2:J-1
    rhonew(j)  = rhoold(j)   - lambda*(roeflux1(j) - roeflux1(j-1));
    mnew(j)    = mold(j)     - lambda*(roeflux2(j) - roeflux2(j-1));
    totenew(j) = totenold(j) - lambda*(roeflux3(j) - roeflux3(j-1));
  end 
  uold   = mnew./rhonew;  	     press = gam1*(totenew - 0.5*mnew.*uold);
  eold   = gammab*press./rhonew;      htot = (totenew + press)./rhonew;
  rhoold = rhonew;  mold = mnew;  totenold = totenew;
  cnew   = sqrt(gamma*gam1*eold);     mach = uold./cnew;
end
entropy = log(press./rhoold.^gamma);

figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),  plot(xcenter,rhonew,'o')
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14), plot(xcenter,uold,'o')
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14), plot(xcenter,press,'o')
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14), plot(xcenter,mach,'o')
subplot(2,3,5),hold on,title('ENTROPY','fontsize',14), plot(xcenter,entropy,'o')
subplot(2,3,6),hold on,title('Roe scheme','fontsize',14), axis('off')

Riemann		% Plot exact solution
