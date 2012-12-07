function [x,w] = Lobatto(n)

%*****************************************************************************80
%
%% LOBATTO_COMPUTE computes a Lobatto quadrature rule.
%
%  Discussion:
%
%    The integral:
%
%      Integral ( -1 <= X <= 1 ) F(X) dX
%
%    The quadrature rule:
%
%      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
%
%    The quadrature rule will integrate exactly all polynomials up to
%    X**(2*N-3).
%
%    The Lobatto rule is distinguished by the fact that both endpoints
%    (-1 and 1) are always abscissas of the rule.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 February 2007
%
%  Author:
%
%    Original MATLAB code by Greg von Winckel.
%    This MATLAB version by John Burkardt.
%
%  Reference:
%
%    Milton Abramowitz, Irene Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964,
%    ISBN: 0-486-61272-4,
%    LC: QA47.A34.
%
%    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
%    Spectral Methods in Fluid Dynamics,
%    Springer, 1993,
%    ISNB13: 978-3540522058,
%    LC: QA377.S676.
%
%    Arthur Stroud, Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966,
%    LC: QA299.4G3S7.
%
%    Daniel Zwillinger, editor,
%    CRC Standard Mathematical Tables and Formulae,
%    30th Edition,
%    CRC Press, 1996,
%    ISBN: 0-8493-2479-3.
%
%  Parameters:
%
%    Input, integer N, the order.  
%    N must be at least 2.
%
%    Output, real X(N), the abscissas.
%
%    Output, real W(N), the weights,
%
 if ( n < 2 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'LOBATTO_COMPUTE - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of N = %d\n', n );
    fprintf ( 1, '  N must be at least 2.\n' );
    error ( 'LOBATTO_COMPUTE - Fatal error!' );
  end

  x = zeros ( n, 1 );
  w = zeros ( n, 1 );

  tolerance = 100.0 * eps;
%
%  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
%
  x(1:n,1) = cos ( pi * ( 0 : n - 1 ) / ( n - 1 ) )';
  xold(1:n,1) = 2.0;

  while ( tolerance < max ( abs ( x(1:n,1) - xold(1:n,1) ) ) )

    xold(1:n,1) = x(1:n,1);

    p(1:n,1) = 1.0;
    p(1:n,2) = x(1:n,1);

    for j = 2 : n-1
      p(1:n,j+1) = ( ( 2 * j - 1 ) * x(1:n,1) .* p(1:n,j)     ...
                   + (   - j + 1 ) *             p(1:n,j-1) ) ...
                   / (     j     );
    end

    x(1:n,1) = xold(1:n,1) - ( x(1:n,1) .* p(1:n,n) - p(1:n,n-1) ) ...
             ./ ( n * p(1:n,n) );
  end

  x(1:n) = r8vec_reverse ( n, x );

  w(1:n,1) = 2.0 ./ ( ( n - 1 ) * n * p(1:n,n).^2 );

  return
end
