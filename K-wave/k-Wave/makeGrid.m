function kgrid = makeGrid(varargin)
%MAKEGRID Create k-Wave grid structure.
%
% DESCRIPTION:
%       makeGrid creates an object of the kWaveGrid class containing the
%       grid coordinates and wavenumber matrices used in k-Wave
%       simulations. The grid structures are indexed as: (x, 1) in 1D; 
%       (x, y) in 2D; and (x, y, z) in 3D.   
%
%       Note, for older versions of MATLAB in which custom class
%       definitions are not supported, makeGrid returns a MATLAB grid
%       structure with the same fields. This structure requires additional
%       memory compared to an object of the kWaveGrid class.   
%
% USAGE:
%       kgrid = makeGrid(Nx, dx)
%       kgrid = makeGrid(Nx, dx, Ny, dy)
%       kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz)
%
% INPUTS:
%       Nx, Ny, Nz  - number of grid points in each Cartesian direction
%       dx, dy, dz  - grid point spacing in each Cartesian direction [m]
%
% OUTPUTS:
%       kgrid       - k-Wave grid structure used by the simulation and
%                     reconstructions functions within k-Wave
%
%       This has the following properties:
%
%       kgrid.k       - ND grid of the scalar wavenumber
%       kgrid.k_max   - maximum spatial frequency supported by the grid
%       kgrid.t_array - evenly spaced array of time values (set to 'auto')
%       kgrid.Nt      - number of time steps (set to 'auto')
%       kgrid.dt      - time step [s] (set to 'auto')
%       kgrid.dim     - number of spatial dimensions
%       kgrid.total_grid_points - total number of grid points
%
%       And for each spatial dimension x, y, z:
%
%       kgrid.Nx      - number of grid points
%       kgrid.dx      - grid point spacing [m]
%       kgrid.x       - plaid ND grid of the x coordinates centered about 0
%                       [m]  
%       kgrid.x_vec   - 1D vector of the x coordinate [m]
%       kgrid.x_size  - length of grid dimension [m]
%       kgrid.kx      - plaid ND grid of the wavenumber components centered
%                       about 0
%       kgrid.kx_vec  - 1D vector of the wavenumber components
%       kgrid.kx_max  - maximum spatial frequency supported by the grid in
%                       the x direction 
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 12th March 2009
%       last update - 21st July 2011
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also cart2grid, interpCartData, kspaceFirstOrder1D,
% kspaceFirstOrder2D, kspaceFirstOrder3D, kspaceSecondOrder, ndgrid,
% makeTime, smooth

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

try
    % create a new instance of the kWaveGrid class
    kgrid = kWaveGrid(varargin{:});
catch %#ok<CTCH>
    % create a grid structure instead
    disp('User defined classes not supported, returning grid structure...');
    kgrid = makeGridStruct(varargin{:});
end
