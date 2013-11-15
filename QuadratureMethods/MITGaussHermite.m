function int = MITGaussHermite(n)
% returns int struct for gaussHermite quadrature with n points
% int.x = points
% int.w = weights
% int.n = number of points
% 
% recall that sum(f(x).*w) =~ int f(x) exp(-x^2) dx from -infty to +infty
%
% see http://mathworld.wolfram.com/Hermite-GaussQuadrature.html
% for algorithm

  Hn = hermitePoly(n);
  Hn_1 = hermitePoly(n-1);
  int.n = n;
  int.x = zeros(n,1); % just to make sure column vectors
  int.w = int.x; 
  int.x = roots(Hn);
  int.w = 2^(n-1)*factorial(n)*sqrt(pi) ./ ...
          (n*polyval(Hn_1,int.x)).^2;
  % sort x
  [int.x order]=sort(int.x);
  int.w=int.w(order); % make w same order
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hk = hermitePoly(n)
% HermitePoly.m by David Terr, Raytheon, 5-10-04
% Given nonnegative integer n, compute the 
% Hermite polynomial H_n. Return the result as a vector whose mth
% element is the coefficient of x^(n+1-m).
% polyval(HermitePoly(n),x) evaluates H_n(x).

  if n==0 
    hk = 1;
  elseif n==1
    hk = [2 0];
  else
    hkm2 = zeros(1,n+1);
    hkm2(n+1) = 1;
    hkm1 = zeros(1,n+1);
    hkm1(n) = 2;

    for k=2:n
      hk = zeros(1,n+1);
      for e=n-k+1:2:n
        hk(e) = 2*(hkm1(e+1) - (k-1)*hkm2(e));
      end
      hk(n+1) = -2*(k-1)*hkm2(n+1);
      if k<n
        hkm2 = hkm1;
        hkm1 = hk;
      end
    end
  end
end % function hermitePoly()