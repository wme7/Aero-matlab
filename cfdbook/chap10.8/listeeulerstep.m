% Liou-Steffen scheme with limiting: Eulerstep in Runge-Kutta time stepping
% Called by Liou_Steffen_MUSCL_scheme

% Theory in Sections 10.3 and 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Function called: extrap  
   
  [U1L, U1R] = extrap(rhostar,limtype);   	% Left and right 
  [U2L, U2R] = extrap(mstar,limtype);		%extrapolated
  [U3L, U3R] = extrap(totenstar,limtype);	%states
  uL = U2L./U1L; uR = U2R./U1R;
  cL = sqrt(gamma*gam1*(U3L./U1L - 0.5*uL.^2));
  cR = sqrt(gamma*gam1*(U3R./U1R - 0.5*uR.^2));
  pL = U1L.*cL.^2/gamma; pR = U1R.*cR.^2/gamma; 
  machL = uL./cL;  machR = uR./cR;
  rhocL = U1L.*cL; rhocR = U1R.*cR;
  enthalpyL = cL.^2/(gamma-1) + 0.5*uL.^2;
  enthalpyR = cR.^2/(gamma-1) + 0.5*uR.^2;
  
% Preallocations
  machplus = machL; machminus = machL; machhalf = machL;
  presplus = machL; presminus = machL;
  
  for j = 1:J-1
    if abs(machL(j)) > 1
      machplus(j) = 0.5*(machL(j) + abs(machL(j)));
      presplus(j) = 0.5*pL(j)*(1 + abs(machL(j))/machL(j));
    else
      machplus(j) = 0.25*(machL(j) + 1)^2;
      presplus(j) = 0.5*pL(j)*(1 + machL(j));
    end
    if abs(machR(j)) > 1
      machminus(j) = 0.5*(machR(j) - abs(machR(j)));
      presminus(j) = 0.5*pR(j)*(1 - abs(machR(j))/machR(j));
    else
      machminus(j) = - 0.25*(machR(j) - 1)^2;
      presminus(j) = 0.5*pR(j)*(1 - machR(j));
    end   
  end
  
% Liou-Steffen fluxes
  machhalf = machplus + machminus;
  flux1a = 0.5*(machhalf + abs(machhalf)).*rhocL;
  flux1b = 0.5*(machhalf - abs(machhalf)).*rhocR;
  flux1  = flux1a + flux1b;
  flux2  = flux1a.*uL + flux1b.*uR + presplus + presminus;
  flux3  = flux1a.*enthalpyL + flux1b.*enthalpyR;
     
% Update of state variables
  rhostar(1) = rholeft; rhostar(J) = rhoright;
  mstar(1) = rholeft*uleft; mstar(J) = rhoright*uright;
  totenstar(1) = totenleft; totenstar(J) = totenright;
  for j = 2:J-1
    rhostar(j)   = rhoold(j)   - rkalpha*lambda*(flux1(j) - flux1(j-1));
    mstar(j)     = mold(j)     - rkalpha*lambda*(flux2(j) - flux2(j-1));
    totenstar(j) = totenold(j) - rkalpha*lambda*(flux3(j) - flux3(j-1));
  end
  
