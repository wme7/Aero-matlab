function [Xl] = advanceBoundary(Xl0, u, X, Y, dt, N, h);
%
% [Xl] = advanceBoundary(Xl0, u, X, Y, dt, N, h);
%
%  Advances one time step the Lagrangian Points Xl0
%
%  Returns:
%     Xl  = the new position of the Lagrangian Points Xl0 after dt
%
%  Input:
%     Xl0 = the old position of the Lagrangian Points
%     u   = the current velocity field
%     X   = N*N length vector of X location of points
%     Y   = ", but Y location of points
%     dt  = time step
%     N   = number of points in each spatial direction
%     h   = mesh width
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


[idxs delta] = evalDelta(Xl0, X, Y, N, h);

u1 = u(:,1);
u2 = u(:,2);
Xl = Xl0 + (dt*h*h) * [ dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2) ];

