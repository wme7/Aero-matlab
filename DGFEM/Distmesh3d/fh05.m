function h = fh05 ( p, varargin )

%*****************************************************************************80
%
%% FH05 returns a mesh size function for the cylinder with a hole.
%
%  Copyright:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Reference:
%
%    Per-Olof Persson and Gilbert Strang,
%    A Simple Mesh Generator in MATLAB,
%    SIAM Review,
%    Volume 46, Number 2, June 2004, pages 329-345.
%
%  Parameters:
%
%    Input, real P(NP,3), the point coordinates.
%
%    Input, VARARGIN, room for extra arguments.
%
%    Output, real H(NP,1), the mesh size function.
%
  h1 = 4.0 * sqrt ( sum ( p.^2, 2 ) ) - 1.0;

  h = min ( h1, 2.0 );

  return
end