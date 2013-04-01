function int = MITGaussLaguerre( norder, alpha )

%% LAGUERRE_COM computes the abscissa and weights for Gauss-Laguerre
%quadrature.
%
%  Discussion:
%
%    In the simplest case, ALPHA is 0, and we are approximating the
%    integral from 0 to INFINITY of EXP(-X) * F(X).  When this is so,
%    it is easy to modify the rule to approximate the integral from
%    A to INFINITY as well.
%
%    If ALPHA is nonzero, then there is no simple way to extend the
%    rule to approximate the integral from A to INFINITY.  The simplest
%    procedures would be to approximate the integral from 0 to A.
%
%    The integration interval is [ A, +Infinity ) or [ 0, +Infinity ).
%
%    The weight function w(x) = EXP ( - X ) or EXP ( - X ) * X**ALPHA.
%
%    The integral to approximate:
%
%      Integral ( A <= X < +INFINITY ) EXP ( - X ) * F(X) dX
%    or
%      Integral ( 0 <= X < +INFINITY ) EXP ( - X ) * X**ALPHA * F(X) dX
%
%    The quadrature rule:
%
%      EXP ( - A ) * Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( A+XTAB(I) )
%    or
%      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
%
%  Modified:
%
%    12 October 2005
%
%  Reference:
%
%    Arthur Stroud and Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966.
%
%  Parameters:
%
%    Input, integer NORDER, the order of the quadrature rule to be
%    computed.
%    NORDER must be at least 1.
%
%    Input, real ALPHA, the exponent of the X factor.
%    Set ALPHA = 0.0 for the simplest rule.
%    ALPHA must be nonnegative.
%
%    Output, real XTAB(NORDER), the Gauss-Laguerre abscissas.
%
%    Output, real WEIGHT(NORDER), the Gauss-Laguerre weights.
%

%
%  Set the recursion coefficients.
%
  for i = 1 : norder
    b(i) = ( alpha + 2 * i - 1 );
  end

  for i = 1 : norder
    c(i) = ( i - 1 ) * ( alpha + i - 1 );
  end

  cc = gamma ( alpha + 1.0 ) * prod ( c(2:norder) );
  int.x = zeros(norder,1);
  int.w = int.x;
  int.n =norder;
  for i = 1 : norder
%
%  Compute an estimate for the root.
%
    if ( i == 1 )

      x = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) ...
        / ( 1.0 + 2.4 * norder + 1.8 * alpha );

    elseif ( i == 2 )

      x = x + ( 15.0 + 6.25 * alpha ) ...
        / ( 1.0 + 0.9 * alpha + 2.5 * norder );

    else

      r1 = ( 1.0 + 2.55 * ( i - 2 ) ) / ( 1.9 * ( i - 2 ) );

      r2 = 1.26 * ( i - 2 ) * alpha / ( 1.0 + 3.5 * ( i - 2 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      x = x + ratio * ( x - xtab(i-2) );

    end
%
%  Use iteration to find the root.
%
    [ x, dp2, p1 ] = laguerre_root ( x, norder, alpha, b, c );
%
%  Set the abscissa and weight.
%
    int.x(i) = x;
    int.w(i) = ( cc / dp2 ) / p1;

  end
  int.n = norder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ x, dp2, p1 ] = laguerre_root ( x, norder, alpha, b, c )

%% LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
%
%  Modified:
%
%    12 October 2005
%
%  Reference:
%
%    Arthur Stroud and Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966.
%
%  Parameters:
%
%    Input, real X, the approximate root, which
%    should be improved on output.
%
%    Input, integer NORDER, the order of the polynomial to be computed.
%
%    Input, real ALPHA, the exponent of the X factor.
%
%    Input, real B(NORDER), C(NORDER), the recursion coefficients.
%
%    Output, real X, the approximate root, which
%    should be improved on output.
%
%    Output, real DP2, the value of L'(NORDER)(X).
%
%    Output, real P1, the value of L(NORDER-1)(X).
%
  maxstep = 10;

  eps = d_epsilon ( x );

  for i = 1 : maxstep

    [ p2, dp2, p1 ] = laguerre_recur ( x, norder, alpha, b, c );

    d = p2 / dp2;
    x = x - d;

    if ( abs ( d ) <= eps * ( abs ( x ) + 1.0 ) )
      return;
    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ p2, dp2, p1 ] = laguerre_recur ( x, norder, alpha, b, c )

%% LAGUERRE_RECUR finds the value and derivative of a Laguerre
%polynomial.
%
%  Modified:
%
%    12 October 2005
%
%  Reference:
%
%    Arthur Stroud and Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966.
%
%  Parameters:
%
%    Input, real X, the point at which polynomials are evaluated.
%
%    Input, integer NORDER, the order of the polynomial to be computed.
%
%    Input, real ALPHA, the exponent of the X factor in the
%    integrand.
%
%    Input, real B(NORDER), C(NORDER), the recursion
%    coefficients.
%
%    Output, real P2, the value of L(NORDER)(X).
%
%    Output, real DP2, the value of L'(NORDER)(X).
%
%    Output, real P1, the value of L(NORDER-1)(X).
%
  p1 = 1.0;
  dp1 = 0.0;

  p2 = x - alpha - 1.0;
  dp2 = 1.0;

  for i = 2 : norder

    p0 = p1;
    dp0 = dp1;

    p1 = p2;
    dp1 = dp2;

    p2 = ( x - b(i) ) * p1 - c(i) * p0;
    dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0;

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%