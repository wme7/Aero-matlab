% Osher scheme (O-variant) for one-dimensional Euler equations

% Copyright 2001 P. Wesseling
% This program and its subprograms may be freely used, modified and distributed
% under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html

% Theory in Section 10.4 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% This program makes Figs. 10.11 -- 10.14 in the book

% Functions called: f, problem_specification, Riemann 

global  PRL  CRL MACHLEFT  gamma  pleft  pright  rholeft  rhoright  uleft...
	uright  tend  lambda		% lambda = dt/dx

		% .....................Input............................
gamma = 1.4; 	% Ratio of specific heats
J = 48;		% Number of grid cells
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
  else, 	     press(j) = pright;  rhoold(j) = rhoright; uold(j) = uright;
  end 
end

	% Initialization of cell center variables
totenold = rhoold.*(0.5*uold.*uold + gammab*press./rhoold); % Total energy rho*E
totenleft = totenold(1); totenright = totenold(J);
mold = rhoold.*uold;					    % Momentum m
c = sqrt(gamma*press./rhoold);				    % Sound speed

	% Preallocation of auxiliary variables
aalpha = zeros(J-1,1); c13 = aalpha; uH = aalpha; p13 = aalpha;

count = 0;			% Counter of flux evaluations
flopcount = flops;		% Counter of floating point operations

t = 0;
for i = 1:n,  t = t + dt;
  lambda1 = uold + c; lambda3 = uold - c;		% Two eigenvalues
  Psi0 = uold - 2*gammab*c; Psi1 = uold + 2*gammab*c;	% Riemann invariants
  z = log(press./(rhoold.^gamma));
  for j = 1:J-1
    aalpha(j) = exp((z(j+1)-z(j))/(2*gamma));
    c13(j) = -0.5*gam1*(Psi0(j)-Psi1(j+1))/(1+aalpha(j));	% c13 = c_{1/3}
    uH(j) = (aalpha(j)*Psi0(j)+Psi1(j+1))/(aalpha(j)+1);
    p13(j) = (gamgam*exp(z(j))*(c13(j)^(-2*gamma)))^(-gammab);
  end    
  lam113 = uH + c13;			% lam113 = eigenvalue_1{_1/3}
  c23 = aalpha.*c13;			% c23 = c_{2/3}  
  lam323 = uH - c23; 			% lam323 = eigenvalue_3{2/3}   
  flux1 = zeros(J-1,1);			% Preallocation of Osher fluxes 
  flux2 = flux1; flux3 = flux1;

  for j = 1:J-1				% Computation of Osher fluxes
    ss = 1 + sign(lambda1(j));
    if ss ~= 0
      count = count + 1;		% Update of flux evaluation counter
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*mold(j);
      flux2(j) = flux2(j) + ssh*(mold(j)^2/rhoold(j) + press(j));
      flux3(j) = flux3(j) +...
        ssh*(uold(j)*(gamma*totenold(j) - 0.5*gam1*mold(j)*uold(j)));
    end
  end
  for j = 1:J-1
    ss = sign(lam113(j)) - sign(lambda1(j));
    if ss ~= 0
      count = count + 1;		% Update of flux evaluation counter
      usonic = (gam1/(gamma+1))*Psi0(j); csonic = - usonic;
      psonic = (gamgam*exp(z(j))/csonic^(2*gamma))^(-gammab);
      rhosonic = gamma*psonic/(csonic^2);
      hsonic = (0.5+gammab)*csonic^2; msonic = rhosonic*usonic;    
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*msonic;
      flux2(j) = flux2(j) + ssh*(psonic + usonic*msonic);
      flux3(j) = flux3(j) + ssh*hsonic*msonic;
    end
  end
  for j = 1:J-1
    ss = sign(uH(j)) - sign(lam113(j));
    if ss ~= 0
      count = count + 1;		% Update of flux evaluation counter
      rho13 = gamma*p13(j)/(c13(j)^2);
      h13 = gammab*(c13(j)^2) + 0.5*(uH(j)^2); m13 = rho13*uH(j);
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*m13;
      flux2(j) = flux2(j) + ssh*(p13(j) + uH(j)*m13);
      flux3(j) = flux3(j) + ssh*h13*m13;
    end
  end
  for j = 1:J-1
    ss = sign(lam323(j)) - sign(uH(j));
    if ss ~= 0
      count = count + 1;		% Update of flux evaluation counter 
      rho23 = gamma*p13(j)/(c23(j)^2);		% NB: p_{2/3} =  p_{1/3}
      h23   = gammab*c23(j)^2 + 0.5*uH(j)^2;   m23 = rho23*uH(j);
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*m23;
      flux2(j) = flux2(j) + ssh*(p13(j) + uH(j)*m23);
      flux3(j) = flux3(j) + ssh*h23*m23;
    end
  end
  for j = 1:J-1
    ss = sign(lambda3(j+1)) - sign(lam323(j));
    if ss~= 0
      count = count + 1;		% Update of flux evaluation counter
      usonic   = (gam1/(gamma+1))*Psi1(j+1);   csonic =  usonic;
      psonic   = (gamgam*exp(z(j+1))*csonic^(-2*gamma))^(-gammab);
      rhosonic = gamma*psonic/(csonic^2);
      hsonic   = (0.5+gammab)*(csonic^2); msonic = rhosonic*usonic;    
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*msonic;
      flux2(j) = flux2(j) + ssh*(psonic + usonic*msonic);
      flux3(j) = flux3(j) + ssh*hsonic*msonic;   
    end
  end
  for j = 1:J-1
    ss = 1 - sign(lambda3(j+1));
    if ss ~= 0,   ssh = 0.5*ss;
      count = count + 1;		% Update of flux evaluation counter
      flux1(j) = flux1(j) + ssh*mold(j+1);
      flux2(j) = flux2(j) + ssh*(mold(j+1)^2/rhoold(j+1) + press(j+1));
      flux3(j) = flux3(j) +...
        ssh*(uold(j+1)*(gamma*totenold(j+1)-0.5*gam1*mold(j+1)*uold(j+1)));
    end
  end
  
	% Update of state variables
  rhonew(1)  = rholeft;		rhonew(J)  = rhoright;
  mnew(1)    = rholeft*uleft;	mnew(J)    = rhoright*uright;
  totenew(1) = totenleft; 	totenew(J) = totenright;
  for j = 2:J-1
    rhonew(j)   = rhoold(j)   - lambda*(flux1(j) - flux1(j-1));
    mnew(j)     = mold(j)     - lambda*(flux2(j) - flux2(j-1));
    totenew(j)  = totenold(j) - lambda*(flux3(j) - flux3(j-1));
  end
  uold   = mnew./rhonew;   	press    = gam1*(totenew - 0.5*mnew.*uold);
  rhoold = rhonew; 		totenold = totenew; mold = mnew;
  c      = sqrt(gamma*press./rhoold); 
end

flopcount = flops - flopcount
number_of_flux_exaluations = count

mach = uold./c;  entropy = log(press./rhoold.^gamma);
figure(1), clf
subplot(2,3,1),hold on,title('DENSITY','fontsize',14),plot(xcenter,rhonew,'o');
subplot(2,3,2),hold on,title('VELOCITY','fontsize',14),plot(xcenter,uold,'o');
subplot(2,3,3),hold on,title('PRESSURE','fontsize',14),plot(xcenter,press,'o');
subplot(2,3,4),hold on,title('MACHNUMBER','fontsize',14),plot(xcenter,mach,'o');
subplot(2,3,5),hold on,title('ENTHALPY','fontsize',14),
enthalpy = (gamma/(gamma-1))*press./rhonew; plot(xcenter,enthalpy,'o');
subplot(2,3,6),hold on,title('ENTROPY','fontsize',14),plot(xcenter,entropy,'o');

Riemann_Osher		% Plot exact solution

