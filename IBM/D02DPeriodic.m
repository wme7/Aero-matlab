function [D0x, D0y] = D02DPeriodic(N, h);
%
%  [D0x, D0y] = D02DPeriodic(N, h);
%
%  Computes the 2D Periodic centered difference operators 
%     on a square.
%
%  Returns:
%     D0x = centered diff. op in the X direction
%     D0y = centered diff. op in the Y direction
%     
%  Input:
%     N = the number of mesh points in each direction
%     h = the mesh width
%
%
%
%
%
%  License: This code is free to use for any purposes, provided
%           any publications resulting from the use of this code
%           reference the original code/author.
%
%  Author:  Samuel Isaacson (isaacson@math.utah.edu)
%  Date:    11/2007
%
%  Please notify the author of any bugs, and contribute any
%  modifications or bug fixes back to the original author.
%
%  Disclaimer:
%   This code is provided as is. The author takes no responsibility 
%   for its results or effects.


M = N * N;
e              = ones(M,1);
eu1            = ones(M,1);
eu1((N+1):N:M) = 0;
ed1            = ones(M,1);
ed1(N:N:M)     = 0;

% Periodic 2D Centered Difference Operators on Square:
D0x  = spdiags([-ed1 eu1], [-1 1], M, M);
D0y  = spdiags([-e e], [-N N], M, M);

for( i = 1:N )
  D0x( (i-1)*N + 1, i*N ) = -1;
  D0x( i*N, (i-1)*N + 1 ) = 1;
  D0y( i, N*N-N+i )       = -1;
  D0y( N*N-N+i, i )       = 1;
end

D0x = D0x / (2*h);
D0y = D0y / (2*h);
