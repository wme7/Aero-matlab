function [f] = spreadForce(F, Xl, X, Y, dtheta, N, h);
%
% [f] = spreadForce(F, Xl, X, Y, N, h);
%
%  Spreads a force on the Lagrangian points, Xl to the fluid 
%     at pts (X,Y)
%  
%  Returns:
%     f = force on fluid 
%
%  Input:
%     F      = force at Lagrangian points
%     Xl     = location of Lagrangian points
%     X      = N*N vector of X location of fluid mesh point
%     Y      = ", but Y locations
%     dtheta = lagrangian point spacing
%     N      = number of mesh points in each direction
%     h      = mesh width
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


[idxs delta] = evalDelta(Xl, X, Y, N, h);

f = zeros(N*N,2);

% for each Langrangian point calculate the force to apply:
for( l = 1:length(Xl(:,1)) )
  f(idxs(l,:),1) = f(idxs(l,:),1) + F(l,1) * dtheta * delta(l,:)';
  f(idxs(l,:),2) = f(idxs(l,:),2) + F(l,2) * dtheta * delta(l,:)';
end

