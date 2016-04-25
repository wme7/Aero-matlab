function [phi] = evalPhi( r );
%
% [phi] = evalPhi( r );
%
%  Calculates phi( r ), phi the discrete delta function (without scaling)
% 
%  Returns:
%     phi = value of unscaled discrete delta function at each point, r
%
%  Input:
%     r   = the points to evaluate the discrete delta function at
%
%  Note: 
%     This code doesn't handle wrapping across periodic boundaries.
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


phi  = zeros(length(r), 1);
aR   = abs( r );
idx1 = find( aR < 1 );
idx2 = find( aR >= 1 );

rT        = aR(idx1);
phi(idx1) = 3 - 2*rT + sqrt( 1 + 4*rT - 4*rT.*rT );

rT        = aR(idx2);
phi(idx2) = 5 - 2*rT - sqrt( -7 + 12*rT - 4*rT.*rT );

phi = phi / 8;

