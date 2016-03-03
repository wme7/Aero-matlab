function pml = getPML(Nx, dx, dt, c, pml_size, pml_alpha, staggered, dimension)
%GETPML   Return the PML variable.
%
% DESCRIPTION:
%       getPML returns a 1D perfectly matched layer variable based on the
%       given size and absorption coefficient.
%
% USAGE:
%       pml = getPML(Nx, dx, dt, c, pml_size, pml_alpha, staggered, dimension)
%
% INPUTS:
%       Nx          - grid size [grid points]
%       dx          - grid spacing [m]
%       dt          - time spacing [s]
%       c           - sound speed [m/s]
%       pml_size    - size of the PML on each side of the grid [grid points]
%       pml_alpha   - absorption coefficient [Nepers per grid point]
%       staggered   - boolean controlling whether the grid is staggered
%       dimension   - direction of the pml vector (1, 2, or 3)
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 8th February 2012
%       last update - 13th February 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

x = 1:pml_size;
if staggered
    % calculate the curved components of the pml using a staggered grid
    pml_left = pml_alpha*(c/dx)* ( ((x + 0.5) - pml_size - 1) ./ (0 - pml_size) ).^4; 
    pml_left = exp(-pml_left*dt/2);
    pml_right = pml_alpha*(c/dx)* ( (x + 0.5) ./ pml_size ).^4;
    pml_right = exp(-pml_right*dt/2);
else
    % calculate the curved components of the pml using a regular grid
    pml_left = pml_alpha*(c/dx)* ( (x - pml_size - 1) ./ (0 - pml_size) ).^4;
    pml_left = exp(-pml_left*dt/2);
    pml_right = pml_alpha*(c/dx)*( x./pml_size ).^4;
    pml_right = exp(-pml_right*dt/2);
end

% add the components of the pml to the total function
pml = ones(1, Nx);
pml(1:pml_size) = pml_left;
pml(end-pml_size+1:end) = pml_right;

% reshape the pml vector to be in the desired direction
switch dimension
    case 1
        pml = pml.';
    case 3
        pml = reshape(pml, [1 1 Nx]);
end

% Other forms:
% ------------
% Use this to include an extra unity point:
% pml_left = pml_alpha*(c/dx)* ( (x - pml_size) ./ (1 - pml_size) ).^2;
% pml_right = pml_alpha*(c/dx)* ( (x - 1) ./ (pml_size - 1) ).^2;
% Staggered grid equivalents:
% pml_left = pml_alpha*(c/dx)* ( ((x + 0.5) - pml_size) ./ (1 - pml_size) ).^2; 
% pml_right = pml_alpha*(c/dx)* ( ((x + 0.5) - 1) ./ (pml_size - 1) ).^2;