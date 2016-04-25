function [L] = lap2DPeriodic(N, h);
%
% [L] = lap2DPeriodic(N, h)
%
%  Constructs the 2D Laplacian for a periodic square mesh using
%     the standard 5-point stencil.
%
%  Returns:
%     L = discrete Laplacian
%
%  Input:
%     N = number of mesh points in each direction
%     h = mesh width
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


% Periodic 2D Laplace Operator on Square:
isPeriodic     = 1;
e              = ones(M,1);
ed             = -4*e;
if( ~isPeriodic )
  ed(N:N:M)      = -3;
  ed(1:N:M)      = -3;
  ed(1)          = -2;
  ed(N)          = -2;
  ed((N-1)*N+1)  = -2;
  ed(N*N)        = -2;
end
eu1            = ones(M,1);
eu1((N+1):N:M) = 0;
ed1            = ones(M,1);
ed1(N:N:M)     = 0;
  
L = spdiags([e ed1 ed eu1 e], [-N -1:1 N], M, M);

if( isPeriodic )
  for(i = 1:N)
    L( i, N*N-N+i )     = 1;
    L( N*N-N+i, i )     = 1;
    L( (i-1)*N+1, i*N ) = 1;
    L( i*N, (i-1)*N+1 ) = 1;
  end
end

L = L ./ (h*h);
