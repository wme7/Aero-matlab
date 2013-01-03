% Staggered unified compressible/incompressible scheme for one-dimensional 
%	Euler equations
% Eulerstep in Runge-Kutta time stepping
% Called by staggered_unified_scheme

% Theory in Section 14.6 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Function called: extrap

allam = rkalpha*lambda; 
rhoL = extrap(rhostar, limtype);
rhostar(1) = rhoold(1) - allam*(ustar(2)*rhoL(1) - uleft*rholeft);
for j = 2:J-1
  rhostar(j) = rhoold(j) - allam*(ustar(j+1)*rhoL(j) - ustar(j)*rhoL(j-1));
end
um = mstar.*ustar;	umL = extrap(um, limtype);
mstar(1) = rholeft*uleft;
for j = 2:J-1
  mstar(j) = mold(j) - allam*(umL(j) - umL(j-1) + pold(j) - pold(j-1));
end
mstar(J) = mold(J) - allam*(um(J) - umL(J-1) + pright - pold(J-1));
rhowall(1) = rholeft; rhowall(J) = rhostar(J-1);
rhoL = extrap(rhostar, limtype);
for j = 2:(J-1)
  rhowall(j) = 0.5*(rhostar(j-1) + rhostar(j));
end
ustar = mstar./rhowall;

