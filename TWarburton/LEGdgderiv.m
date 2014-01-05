function [dfdx] = LEGdgderiv(f, nodex, tauL, tauR, D, F, G, H, J)

  [p,N] = size(f);

  % p = maximum polynomial order 
  p=p-1;

  % N = number of nodes
  N=N+1;

  % space for derivative
  dfdx = zeros(p+1,N-1);

  % cells 2 to N-2  
  ids = (2:N-2);
  dfdx(:,ids) = (D*f(:,ids) + tauL*F*f(:,ids) + tauL*G*f(:,ids-1) + tauR*H*f(:,ids) + tauR*J*f(:,ids+1));

  % cell 1
  dfdx(:,1)   = (D*f(:,1)   + tauL*F*f(:,1)   + tauL*G*f(:,N-1)   + tauR*H*f(:,1)   + tauR*J*f(:,2));

  % cell N-1
  dfdx(:,N-1) = (D*f(:,N-1) + tauL*F*f(:,N-1) + tauL*G*f(:,N-2)   + tauR*H*f(:,N-1) + tauR*J*f(:,1));

  % apply chain rule for physical cell width
  dx    = nodex(2:N)-nodex(1:N-1);
  coeff = ones(p+1,1)*(2./dx);
  dfdx  = dfdx.*coeff;
