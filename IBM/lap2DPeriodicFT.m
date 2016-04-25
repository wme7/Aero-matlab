function [LT] = lap2DPeriodicFT(N, h);
%
% [LT] = lap2DPeriodicFT(N, h);
%  
%  Computes the Fourier Transform of the 5-point 2D periodic Laplacian
%     on a square.
%
%  Returns:
%     LT = 2D FT of the Laplacian as a vector
%
%  Input:
%     N  = number of mesh points in each direction
%     h  = mesh width
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


% 2D Fourier Space Indexes:
l    = repmat( (0:(N-1))', N, 1 );
m    = reshape( repmat( (0:(N-1)), N, 1 ), M, 1);


% Fourier Transformed 2D Laplace Operator on Square:
LT   = -(4 / (h*h)) * ( sin( (pi/N) * l ).^2 + sin( (pi/N) * m ).^2 );
