function s=superbee(a,b)
% The superbee limiter
  s=maxmod(minmod(a,2*b),minmod(2*a,b));