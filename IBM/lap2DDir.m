function [L] = lap2DDir(N, h);
%
% [L] = lap2DDir(N, h)
%
%  Constructs the 2D Laplacian for a square mesh with zero Dirichlet 
%     boundary conditions using the standard 5-point stencil.
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


% 2D Laplace Operator on Square:
e              = ones(M,1);
ed             = -4*e;
eu1            = ones(M,1);
eu1((N+1):N:M) = 0;
ed1            = ones(M,1);
ed1(N:N:M)     = 0;
  
L = spdiags([e ed1 ed eu1 e], [-N -1:1 N], M, M);


L = L ./ (h*h);
