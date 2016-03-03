function [t_array, dt] = makeTime(kgrid, c, cfl, t_end)
%MAKETIME   Create an evenly spaced array of time points.
%
% DESCRIPTION:
%       makeTime creates an evenly spaced array of time points for use with
%       the k-Wave simulation codes based on the Courant-Friedrichs-Lewy
%       stability level cfl and the grid size. The time-step dt is chosen
%       based on the cfl level (the default setting is 0.3), and the number
%       of time-steps is chosen based on the time it takes to travel from
%       one corner of the grid specified by kgrid to the geometrically
%       opposite corner. Note, if c is given as a matrix, the calculation
%       for dt is based on the maximum value, and the calculation for t_end
%       based on the minimum value.
%
% POSSIBLE USAGE:
%       [t_array, dt] = makeTime(kgrid, c)
%       [t_array, dt] = makeTime(kgrid, c, cfl)
%       [t_array, dt] = makeTime(kgrid, c, cfl, t_end)
%       [t_array, dt] = makeTime(kgrid, c, [], t_end)
%
% INPUTS:
%       kgrid       - k-Wave grid structure returned by makeGrid
%       c           - sound speed in the medium [m/s]
%
% OPTIONAL INPUTS:
%       cfl         - Courant-Friedrichs-Lewy stability criterion 
%                     (default = 0.3)
%       t_end       - maximum time [s]
%
% OUTPUTS:
%       t_array     - array of evenly-spaced time points [s]
%       dt          - time step [s]
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 3rd July 2009
%       last update - 2nd September 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also makeGrid, kspaceFirstOrder1D, kspaceFirstOrder2D,
% kspaceFirstOrder3D 

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

% CHANGE LOG:
% 2009/10/05    - changed pixel_dim to select min pixel size instead of max
% 2010/04/09    - added tmax input
% 2011/07/19    - updated with reordered inputs
% 2012/09/02    - modified t_end to use c_min

% check user defined value for the Courant-Friedrichs-Lewy stability
% level, otherwise assign default value
if nargin < 3 || isempty(cfl)
    cfl = 0.3;
end

% if c is a matrix, find the minimum and maximum values
c_max = max(c(:));
c_min = min(c(:));

% check for user define t_end, otherwise set the simulation length based on
% the size of the grid diagonal and the maximum sound speed in the medium
if nargin < 4
    switch kgrid.dim
        case 1
            t_end = kgrid.x_size/c_min;
        case 2
            t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2)/c_min;
        case 3
            t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2 + kgrid.z_size.^2)/c_min;
    end
end

% extract the smallest pixel dimension
switch kgrid.dim
    case 1
        pixel_dim = kgrid.dx;
    case 2
        pixel_dim = min([kgrid.dx, kgrid.dy]);
    case 3
        pixel_dim = min([kgrid.dx, kgrid.dy, kgrid.dz]);
end

% assign a time step based on Courant stability criterion
dt = cfl*pixel_dim/c_max;

% create the time array
t_array = 0:dt:t_end;