function [X W] = hermquad(N)
%
% [X W] = HERMQUAD(N)
%
% Find the Gauss-Hermite abscissae and weights.
%
% Arguments:
%  N - The number of abscissae and weights to return.
%
% Return Values:
%  X - A column vector containing the abscissae.
%  W - A column vector containing the corresponding weights.
%
% Gauss-Hermite quadrature approximates definite integrals of the form
%
%     \int^{-\infty}_{\infty} dx W(x) f(x)
%
% where
%
%     W(x) = \exp( - x^2 )
%
% with the sum
%
%     \sum_{n=1}^{N} w_{n} f(x_{n}).
%
% This function returns the set of abscissae and weights
%
%     {x_{n}, w_{n}}^{N}_{n=1}
%
% for performing this calculation given N, the number of abscissae.
% These abscissae correspond to the zeros of the Nth Hermite
% polynomial.  It can be shown that such integration is exact when f(x)
% is a polynomial of maximum order 2N-1.
%
% The procedure in this calculation is taken more or less directly from
%
% @BOOK{ press-etal-1992a,
%	 AUTHOR   = { Press, William  H.   and 
%                     Flannery, Brian  P.  and
%                     Teukolsky, Saul  A.  and
%                     Vetterling, William  T. },
%       ISBN      = {0521431085},
%	MONTH     = {October},
%	PUBLISHER = {{Cambridge University Press}},
%	TITLE     = {Numerical Recipes in C : The Art of Scientific Computing},
%	YEAR      = {1992}
%      }
%

% precision
EPS = 3.0e-14;

% 1/\pi^{1/4}
PIM4 = 0.7511255444649425;

% maximum number of loops
MAXIT = 10;

% allocate the return values
X = zeros([N 1]);
W = zeros([N 1]);

for i=1:(N+1)/2
  
  % good guesses at initial values for specific roots
  if i == 1
    z = sqrt(2.0*N+1.0) - 1.85575*((2.0*N+1)^(-0.16667));
  elseif i == 2
    z = z - (1.14 * N^0.426 / z);
  elseif i == 3
    z = 1.86 * z - 0.86 * X(1);
  elseif i == 4
    z = 1.91 * z - 0.91 * X(2);
  else
    z = 2.0*z - X(i-2);
  end
  
  for iter=1:MAXIT+1
    p1 = PIM4;
    p2 = 0.0;
    
    for j=1:N
      p3 = p2;
      p2 = p1;
      p1 = z * sqrt(2.0/j) * p2 - sqrt((j-1.0)/j) * p3;
    end
    
    % the derivative
    pp = sqrt(2.0*N) * p2;
    
    % newton step
    z1 = z;
    z  = z1 - p1/pp;
    
    if abs(z-z1) <= EPS
      break;
    end    
  end
  
  if iter == MAXIT+1
    fprintf('Too many iterations in hermquad.\n');
  end
  
  X(i)     = z;
  X(N+1-i) = -z;
  W(i)     = 2.0/(pp*pp);
  W(N+1-i) = W(i);

end
