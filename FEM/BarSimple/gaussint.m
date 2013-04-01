function result = gaussint(P,a,b,ngp)
% Integration using Gauss Quadratures
% by Manuel Diaz 2013.03.20

% compute gauss locations and weights
[w,psi] = gauss1d(ngp);

% compute the change the integration range to -1 to 1
x = (b+a)/2 + (b-a)/2*psi;

% perform quadrature
result  = (b-a)/2*sum(w.*P(x));