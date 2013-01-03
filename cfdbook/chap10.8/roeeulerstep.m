% Roe scheme, Eulerstep in Runge-Kutta time stepping
%   MUSCL for higher order
% Called by Roe_MUSCL_scheme
% Theory in Sections 10.4 and 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Functions called: extrap,  primtoriem,  riemtoprim

global epsi  invariants

  if invariants == 0	% Extrapolation of primitive variables
    [U1L, U1R] = extrap(rhostar, limtype);
    [U2L, U2R] = extrap(mstar, limtype);
    [U3L, U3R] = extrap(totenstar, limtype);
  else			% Extrapolation of Riemann invariants
	% Computation of Riemann invariants
    [Z1, Z2, Z3] = primtoriem(rhostar, mstar, totenstar);
	% Extrapolation of Riemann invariants
    [Z1L, Z1R] = extrap(Z1, limtype); [Z2L, Z2R] = extrap(Z2, limtype);
    [Z3L, Z3R] = extrap(Z3, limtype);
	% Computation of primitive variables
    [U1L, U2L, U3L] = riemtoprim(Z1L, Z2L, Z3L);
    [U1R, U2R, U3R] = riemtoprim(Z1R, Z2R, Z3R);
  end
  
  uL = U2L./U1L; 	uR = U2R./U1R;
  HL =gamma*U3L./U1L - 0.5*(gamma-1)*uL.^2;
  HR =gamma*U3R./U1R - 0.5*(gamma-1)*uR.^2;
     
	% Roe averages
   dd = sqrt(U1R./U1L); hav = (HL + dd.*HR)./(1+dd);
   uav = (uL + dd.*uR)./(1+dd);

   delrho = U1R - U1L;   delm = U2R - U2L;   deltoten = U3R - U3L;

   f1av = 0.5*(U2L + U2R);
   f2av = 0.5*(gam1*U3L - 0.5*(gamma-3)*(U2L.^2)./U1L...
       + gam1*U3R - 0.5*(gamma-3)*(U2R.^2)./U1R);
   f3av = 0.5*(U2L.*HL + U2R.*HR);
   
   cav = sqrt(gam1*(hav - 0.5*uav.*uav));   mav = uav./cav;
   alpha1 = 0.25*mav.*(2+gam1*mav).*delrho -...
           0.5*(1+gam1*mav).*delm./cav + 0.5*gam1*deltoten./(cav.^2);
   alpha2 = (1-0.5*gam1*(mav.^2)).*delrho +...
            gam1*(mav./cav).*delm - gam1*deltoten./(cav.^2);
   alpha3 = -0.25*mav.*(2-gam1*mav).*delrho +...
           0.5*(1-gam1*mav).*delm./cav + 0.5*gam1*deltoten./(cav.^2);
   alambda1 = abs(uav - cav); alambda2 = abs(uav); alambda3 = abs(uav + cav);

	% Harten's sonic entropy fix
  for j = 1:J-1
    if alambda1(j) < epsi
      alambda1(j) = 0.5*(epsi + alambda1(j)^2/epsi);
    end
    if alambda3(j) < epsi
      alambda3(j) = 0.5*(epsi + alambda3(j)^2/epsi);
    end
  end
  
  roeflux1 = f1av - 0.5*alambda1.*alpha1...
     - 0.5*alambda2.*alpha2 - 0.5*alambda3.*alpha3;
  roeflux2 = f2av - 0.5*alambda1.*alpha1.*(uav-cav) -...
     0.5*alambda2.*alpha2.*uav - 0.5*alambda3.*alpha3.*(uav+cav);
  roeflux3 = f3av - 0.5*alambda1.*alpha1.*(hav - uav.*cav) -...
     0.25*alambda2.*alpha2.*uav.*uav - 0.5*alambda3.*alpha3.*(hav + uav.*cav);
  rhostar(1) = rholeft; rhostar(J) = rhoright;
  mstar(1) = rholeft*uleft; mstar(J) = rhoright*uright;
  totenstar(1) = totenleft; totenstar(J) = totenright;
  for j = 2:J-1
    rhostar(j)   = rhoold(j)   - rkalpha*lambda*(roeflux1(j) - roeflux1(j-1));
    mstar(j)     = mold(j)     - rkalpha*lambda*(roeflux2(j) - roeflux2(j-1));
    totenstar(j) = totenold(j) - rkalpha*lambda*(roeflux3(j) - roeflux3(j-1));
  end 
