% Osher scheme, O-variant, Eulerstep in Runge-Kutta time stepping 
%   MUSCL for higher order
% Called by Osher_MUSCL_scheme
% Theory in Sections 10.3 and 10.8 of:

% 	P. Wesseling: Principles of Computational Fluid Dynamics
% 	Springer, Heidelberg, 2000 ISBN 3-5453-0. XII, 642 pp.
% 	See http://ta.twi.tudelft.nl/nw/users/wesseling/cfdbook.html

% Functions called: extrap, primtoriem, riemtoprim  

global invariants

  if invariants == 0	% 1. Extrapolation of primitive variables
    [U1L, U1R] = extrap(rhostar, limtype);
    [U2L, U2R] = extrap(mstar, limtype);
    [U3L, U3R] = extrap(totenstar, limtype);
	% Computation of Riemann invariants
    [Z1L, Z2L, Z3L] = primtoriem(U1L, U2L, U3L);
    [Z1R, Z2R, Z3R] = primtoriem(U1R, U2R, U3R);
  else			% 2. Extrapolation of Riemann invariants
	% Computation of Riemann invariants
    [Z1, Z2, Z3] = primtoriem(rhostar, mstar, totenstar);
	% Extrapolation of Riemann invariants
    [Z1L, Z1R] = extrap(Z1, limtype); [Z2L, Z2R] = extrap(Z2, limtype);
    [Z3L, Z3R] = extrap(Z3, limtype);
	% Computation of primitive variables
    [U1L, U2L, U3L] = riemtoprim(Z1L, Z2L, Z3L);
    [U1R, U2R, U3R] = riemtoprim(Z1R, Z2R, Z3R);
  end
  
  uL = 0.5*(Z1L + Z3L); uR = 0.5*(Z1R + Z3R);
  cL = 0.25*gam1*(Z3L - Z1L); cR = 0.25*gam1*(Z3R - Z1R);  
  pL = U1L.*cL.^2/gamma; pR = U1R.*cR.^2/gamma;
 
  lambda1 = uL+cL; lambda3 = uR-cR;	% Two eigenvalues
  
  aalpha = exp((Z2R-Z2L)/(2*gamma));
  c13  = -0.5*gam1*(Z1L-Z3R)./(1+aalpha);	% c13 = c_{1/3} 
  uH = (aalpha.*Z1L+Z3R)./(aalpha+1);
  p13 = (gamgam*exp(Z2L).*(c13.^(-2*gamma))).^(-gammab);
    
  lam113 = uH + c13;		% lam113 = lambda_1(1/3)
  c23 = aalpha.*c13;		% c23 = c_{2/3}  
  lam323 = uH - c23;		% lam323 = lambda_3(2/3)
      
	% Preallocation of Osher flux
  flux1 = zeros(J-1,1); flux2 = flux1; flux3 = flux1;

  for j = 1:J-1			% Computation of Osher flux
    ss = 1 + sign(lambda1(j));
    if ss ~= 0
      count = count + 1;	% Update of flux counter
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*U2L(j);
      flux2(j) = flux2(j) + ssh*(U2L(j)*uL(j) + pL(j));
      flux3(j) = flux3(j) +...
        ssh*(uL(j)*(gamma*U3L(j) - 0.5*gam1*U2L(j)*uL(j)));
    end
  end
  for j = 1:J-1
   ss = sign(lam113(j)) - sign(lambda1(j));
   if ss ~= 0
      count = count + 1;
      usonic = (gam1/(gamma+1))*Z1L(j); csonic = - usonic;
      psonic = (gamgam*exp(Z2L(j))/csonic^(2*gamma))^(-gammab);
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
      count = count + 1;
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
      count = count + 1; 
      rho23 = gamma*p13(j)/(c23(j)^2);		% p_{2/3} =  p_{1/3}
      h23 = gammab*c23(j)^2 + 0.5*uH(j)^2; m23 = rho23*uH(j);
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*m23;
      flux2(j) = flux2(j) + ssh*(p13(j) + uH(j)*m23);
      flux3(j) = flux3(j) + ssh*h23*m23;
    end
  end
  for j = 1:J-1
    ss = sign(lambda3(j)) - sign(lam323(j));
    if ss~= 0
      count = count + 1;
      usonic = (gam1/(gamma+1))*Z3R(j); csonic =  usonic;
      psonic = (gamgam*exp(Z2R(j))*csonic^(-2*gamma))^(-gammab);
      rhosonic = gamma*psonic/(csonic^2);
      hsonic = (0.5+gammab)*(csonic^2); msonic = rhosonic*usonic;    
      ssh = 0.5*ss;
      flux1(j) = flux1(j) + ssh*msonic;
      flux2(j) = flux2(j) + ssh*(psonic + usonic*msonic);
      flux3(j) = flux3(j) + ssh*hsonic*msonic;   
    end
  end
  for j = 1:J-1
    ss = 1 - sign(lambda3(j));
    if ss ~= 0
      ssh = 0.5*ss;
      count = count + 1;
      flux1(j) = flux1(j) + ssh*U2R(j);
      flux2(j) = flux2(j) + ssh*(U2R(j)*uR(j) + pR(j));
      flux3(j) = flux3(j) +...
        ssh*(uR(j)*(gamma*U3R(j)-0.5*gam1*U2R(j)*uR(j)));
    end
  end
  
	% Update of state variables
  rhostar(1) = rholeft; rhostar(J) = rhoright;
  mstar(1) = rholeft*uleft; mstar(J) = rhoright*uright;
  totenstar(1) = totenleft; totenstar(J) = totenright;
  for j = 2:J-1
    rhostar(j)   = rhoold(j)   - rkalpha*lambda*(flux1(j) - flux1(j-1));
    mstar(j)     = mold(j)     - rkalpha*lambda*(flux2(j) - flux2(j-1));
    totenstar(j) = totenold(j) - rkalpha*lambda*(flux3(j) - flux3(j-1));
  end
  
