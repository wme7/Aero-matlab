function s=mc(a,b)
% The MacCormack limiter.
  h1=0.5*(a+b); 
  h2=2*minmod(a,b);
  s=minmod(h1,h2);