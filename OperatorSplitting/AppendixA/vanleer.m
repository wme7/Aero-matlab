function s=vanleer(a,b)
% The Van-Leer limiter
  small=1e-12;
  aa=sqrt(small+a.^2);
  bb=sqrt(small+b.^2);
  s=b.*a.*(sign(a)+sign(b))./(aa+bb);