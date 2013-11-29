function d = fd04 ( p )

%*****************************************************************************80
%
%% FD04 is a signed distance function for problem 4.
%
%  Modified:
%
%    11 September 2005
%
%  Parameters:
%
%    Input, real P(N,3), one or more points.
%
%    Output, real D(N), the signed distance of each point to the boundary of the region.
%
  d = sqrt ( sum ( p .^ 2, 2 ) ) - 1.0;

  return
end