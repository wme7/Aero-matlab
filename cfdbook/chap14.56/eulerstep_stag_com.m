% Staggered compressible scheme for one-dimensional Euler equations
% Eulerstep in Runge-Kutta time stepping
% Called by staggered_compressible_scheme

% Theory in Section 14.5 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Function called: extrap
  
totenthalpy = etotstar + pstar;  totenthL = extrap(totenthalpy, limtype);
rhoL        = extrap(rhostar, limtype);
etotstar(1) = etotold(1) - rkalpha*lambda*(ustar(2)*totenthL(1)...
              - uleft*(gamma*pleft/(gamma-1) + 0.5*rholeft*uleft^2));
rhostar(1)  = rhoold(1) - rkalpha*lambda*(ustar(2)*rhoL(1) - uleft*rholeft);
for j = 2:J-1
  etotstar(j) = etotold(j) - rkalpha*lambda*(ustar(j+1)*totenthL(j)...
                - ustar(j)*totenthL(j-1));
  rhostar(j)  = rhoold(j) - rkalpha*lambda*(ustar(j+1)*rhoL(j)...
                - ustar(j)*rhoL(j-1));
end
um       = ustar.*mstar;	umL      = extrap(um, limtype);
ustar(1) = uleft;	  	mstar(1) = rholeft*uleft;
for j = 2:J-1
  mstar(j) = mold(j) - rkalpha*lambda*(umL(j) - umL(j-1))...
             - rkalpha*lambda*(pstar(j) - pstar(j-1));
  ustar(j) = 2*mstar(j)/(rhostar(j) + rhostar(j-1));
end
mstar(J) = mold(J) - rkalpha*lambda*(umL(J) - umL(J-1))...
           - rkalpha*lambda*(pright - pstar(J-1));
ustar(J) = 2*mstar(J)/(rhostar(J-1) + rhoright);
for j = 1:J-1
  pstar(j) = (gamma-1)*(etotstar(j) - 0.25*(mstar(j)*ustar(j) +...
             mstar(j+1)*ustar(j+1)));
end


