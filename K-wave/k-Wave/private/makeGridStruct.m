function kgrid = makeGridStruct(varargin)
%MAKEGRID   Create k-Wave grid structure.
%
% DESCRIPTION:
%       makeGrid creates a MATLAB structure containing the grid coordinates
%       and wavenumber matrices for use in k-space simulations and
%       reconstructions. This function is called by makeGrid.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 12th March 2009
%       last update - 4th October 2012
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

% assign the input values to the grid object
kgrid.Nx = varargin{1};
kgrid.dx = varargin{2};
if nargin == 6
    kgrid.Ny = varargin{3};
    kgrid.dy = varargin{4};
    kgrid.Nz = varargin{5};
    kgrid.dz = varargin{6};
elseif nargin == 4
    kgrid.Ny = varargin{3};
    kgrid.dy = varargin{4};    
    kgrid.Nz = 0;
    kgrid.dz = 0;
elseif nargin == 2
    kgrid.Ny = 0;
    kgrid.dy = 0;
    kgrid.Nz = 0;
    kgrid.dz = 0;    
else
    error('Incorrect number of input arguments');
end

% initialise the remaining structure parameters
kgrid.kx_vec = 0;
kgrid.ky_vec = 0;
kgrid.kz_vec = 0;           
kgrid.k = 0;
kgrid.kx_max = 0;
kgrid.ky_max = 0;
kgrid.kz_max = 0;
kgrid.k_max = 0;
kgrid.x = 0;
kgrid.y = 0;
kgrid.z = 0;
kgrid.kx = 0;
kgrid.ky = 0;
kgrid.kz = 0;
kgrid.x_vec = 0;
kgrid.y_vec = 0;
kgrid.z_vec = 0; 
kgrid.x_size = 0;
kgrid.z_size = 0;
kgrid.y_size = 0;
kgrid.total_grid_points = 0;
kgrid.nonuniform = false;

switch nargin
    case 2
        % assign the grid parameters for the x spatial direction
        [nx, kgrid.kx_vec, kgrid.x_size] = makeDim(kgrid.Nx, kgrid.dx);
        
        % define a spatial grid that is centered about 0
        kgrid.x = nx*kgrid.x_size; 
        
        % create the vector position variables
        kgrid.x_vec = kgrid.x;        

        % define wavenumber component centered about 0
        kgrid.kx = kgrid.kx_vec;
        
        % define the scalar wavenumber based on the wavenumber components
        kgrid.k = abs(kgrid.kx);
        
        % define maximum supported frequency
        kgrid.kx_max = max(abs(kgrid.kx(:)));
        kgrid.k_max = kgrid.kx_max;
        
        % set the number of dimensions
        kgrid.dim = 1;
        
        % set the total number of gridpoints
        kgrid.total_grid_points = kgrid.Nx;
    case 4
        % assign the grid parameters for the x and z spatial directions
        [nx, kgrid.kx_vec, kgrid.x_size] = makeDim(kgrid.Nx, kgrid.dx);
        [ny, kgrid.ky_vec, kgrid.y_size] = makeDim(kgrid.Ny, kgrid.dy);
        
        % define a spatial grid that is centered about 0
        [kgrid.x, kgrid.y] = ndgrid(nx*kgrid.x_size, ny*kgrid.y_size); 
        
        % create the vector position variables
        kgrid.x_vec = nx*kgrid.x_size;
        kgrid.y_vec = ny*kgrid.y_size;           

        % define plaid grids of the wavenumber components centered about 0
        [kgrid.kx, kgrid.ky] = ndgrid(kgrid.kx_vec, kgrid.ky_vec);

        % define the scalar wavenumber based on the wavenumber components
        kgrid.k = sqrt(kgrid.kx.^2 + kgrid.ky.^2);
        
        % define maximum supported frequency
        kgrid.kx_max = max(abs(kgrid.kx(:)));
        kgrid.ky_max = max(abs(kgrid.ky(:)));
        kgrid.k_max = min([kgrid.kx_max, kgrid.ky_max]);      
        
        % set the number of dimensions
        kgrid.dim = 2;
        
        % set the total number of gridpoints
        kgrid.total_grid_points = kgrid.Nx*kgrid.Ny;        
    case 6
        % assign the grid parameters for the x ,y and z spatial directions
        [nx, kgrid.kx_vec, kgrid.x_size] = makeDim(kgrid.Nx, kgrid.dx);
        [ny, kgrid.ky_vec, kgrid.y_size] = makeDim(kgrid.Ny, kgrid.dy);        
        [nz, kgrid.kz_vec, kgrid.z_size] = makeDim(kgrid.Nz, kgrid.dz);

        % define a spatial grid that is centered about 0
        [kgrid.x, kgrid.y, kgrid.z] = ndgrid(nx*kgrid.x_size, ny*kgrid.y_size, nz*kgrid.z_size);

        % create the vector position variables
        kgrid.x_vec = nx*kgrid.x_size;
        kgrid.y_vec = ny*kgrid.y_size;          
        kgrid.z_vec = nz*kgrid.z_size;         
        
        % define plaid grids of the wavenumber components centered about 0
        [kgrid.kx, kgrid.ky, kgrid.kz] = ndgrid(kgrid.kx_vec, kgrid.ky_vec, kgrid.kz_vec);

        % define the scalar wavenumber based on the wavenumber components
        kgrid.k = sqrt(kgrid.kx.^2 + kgrid.ky.^2 + kgrid.kz.^2);
        
        % define maximum supported frequency
        kgrid.kx_max = max(abs(kgrid.kx(:)));
        kgrid.ky_max = max(abs(kgrid.ky(:)));
        kgrid.kz_max = max(abs(kgrid.kz(:)));
        kgrid.k_max = min([kgrid.kx_max, kgrid.ky_max, kgrid.kz_max]); 
        
        % set the number of dimensions
        kgrid.dim = 3;      
        
        % set the total number of gridpoints
        kgrid.total_grid_points = kgrid.Nx*kgrid.Ny*kgrid.Nz;        
end

% set t_array to 'auto' by default
kgrid.t_array = 'auto';

% subfunction to create the grid parameters for a single spatial direction
function [nx, kx_vec, x_size] = makeDim(Nx, dx)

% define the discretisation of the spatial dimension such that there is
% always a DC component
if rem(Nx, 2) == 0
    % grid dimension has an even number of points
    nx = ((-Nx/2:Nx/2-1)/Nx).';
else
    % grid dimension has an odd number of points
    nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx).';
end

% force middle value to be zero in case 1/Nx is a recurring number and the
% series doesn't give exactly zero 
nx(floor(Nx/2) + 1) = 0;

% assign the size parameter
x_size = dx*Nx;

% define the wavenumber vector components
kx_vec = (2*pi/dx).*nx;