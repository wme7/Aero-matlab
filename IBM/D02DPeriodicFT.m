function [D0xT, D0yT] = D02DPeriodicFT(N, h);
%
% [D0xT, D0yT] = D02DPeriodicFT(N, h);
%  
%  Computes the Fourier Transform of the 2D periodic centered difference
%     operator on a square.
%
%  Returns:
%     D0xT = 2D FT of D0x, the X direction centered difference op.
%     D0yT = 2D FT of D0y, the Y direction centered difference op
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

% Fourier Transformed 2D Centered Difference Operators on Square:
D0xT = (sqrt(-1) / h) * sin( (2 * pi / N) * l );
D0yT = (sqrt(-1) / h) * sin( (2 * pi / N) * m );