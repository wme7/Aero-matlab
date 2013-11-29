function d = fd05 ( p )

%*****************************************************************************80
%
%% FD05 is a signed distance function for the cylinder with a hole.
%
%  Modified:
%
%    15 September 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real P(N,3), one or more points.
%
%    Output, real D(N), the signed distance of each point to the boundary of the region.
%
  r = sqrt ( p(:,1).^2 + p(:,2).^2 );
  z = p(:,3);

  d1 = r - 1.0;
  d2 = z - 1.0;
  d3 = - z - 1.0;
  d4 = sqrt ( d1.^2 + d2.^2 );
  d5 = sqrt ( d1.^2 + d3.^2 );

  d = dintersect ( dintersect ( d1, d2 ), d3 );
  ix = ( 0.0 < d1 ) & ( 0.0 < d2 );
  d(ix) = d4(ix);
  ix = ( 0.0 < d1 ) & ( 0.0 < d3 );
  d(ix) = d5(ix);

  d = ddiff ( d, dsphere ( p, 0.0, 0.0, 0.0, 0.5 ) );

  return
end