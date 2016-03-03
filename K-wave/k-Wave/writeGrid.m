function writeGrid(filename, grid_size, grid_spacing, pml_size, pml_alpha, Nt, dt, c_ref) 
%WRITEGRID    Write grid and PML properties to a k-Wave HDF5 file.
%
% DESCRIPTION:
%       writeGrid creates and writes the wavenumber grids and PML variables
%       required by the k-Wave C++ code to the HDF5 file specified by the
%       user. 
%
%       List of parameters that are written:
%           Nx
%           Ny
%           Nz
%           Nt
%           dt
%           dx
%           dy
%           dz
%           c_ref
%           ddx_k_shift_pos_r
%           ddx_k_shift_neg_r
%           ddy_k_shift_pos
%           ddy_k_shift_neg
%           ddz_k_shift_pos
%           ddz_k_shift_neg
%           x_shift_neg_r
%           y_shift_neg_r
%           z_shift_neg_r
%           pml_x_sgx
%           pml_y_sgy
%           pml_z_sgz
%           pml_x
%           pml_y
%           pml_z
%           pml_x_alpha
%           pml_y_alpha
%           pml_z_alpha
%           pml_x_size
%           pml_y_size
%           pml_z_size
%
% USAGE:
%       writeGrid(filename, grid_size, grid_spacing, pml_size, pml_alpha, Nt, dt, c_ref) 
%
% INPUTS:
%       filename            - filename and location of the input HDF5 file
%       grid_size           - [Nx, Ny, Nz]
%       grid_spacing        - [dx, dy, dz]
%       pml_size            - [pml_x_size, pml_y_size, pml_z_size]
%       pml_alpha           - [pml_x_alpha, pml_y_alpha, pml_z_alpha]
%       Nt                  - number of time points
%       dt                  - time step
%       c_ref               - scalar sound speed used in the k-space
%                             operator and to define the pml variables
%
% ABOUT:
%       author              - Bradley Treeby
%       date                - 30th May 2013
%       last update         - 21st August 2014
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also h5writeatt, writeAttributes, writeFlags, writeMatrix

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

%#ok<*INUSL>
%#ok<*NASGU>

% get literals
getH5Literals;

% unpack grid size inputs to make code easier to read
Nx          = grid_size(1);
Ny          = grid_size(2);
Nz          = grid_size(3);
dx          = grid_spacing(1);
dy          = grid_spacing(2);
dz          = grid_spacing(3);
pml_x_size  = pml_size(1);
pml_y_size  = pml_size(2);
pml_z_size  = pml_size(3);
pml_x_alpha = pml_alpha(1);
pml_y_alpha = pml_alpha(2);
pml_z_alpha = pml_alpha(3);

% =========================================================================
% CREATE WAVENUMBER AND PML VECTORS
% =========================================================================

% create the wavenumber grids (assuming Nx, Ny and Nz are even)
nx = ((-Nx/2:Nx/2-1)/Nx).';
nx(floor(Nx/2) + 1) = 0;
kx_vec = (2*pi/dx).*nx; 

ny = ((-Ny/2:Ny/2-1)/Ny).';
ny(floor(Ny/2) + 1) = 0;
ky_vec = (2*pi/dy).*ny; 

nz = ((-Nz/2:Nz/2-1)/Nz).';
nz(floor(Nz/2) + 1) = 0;
kz_vec = (2*pi/dz).*nz; 

% force the vector operators be in the correct direction (Nx, 1, 1), (1, Ny, 1), (1, 1, Nz) 
ky_vec = ky_vec.'; 
kz_vec = permute(kz_vec, [2 3 1]);

% create vector derivative and shift variables
ddx_k_shift_pos = ifftshift( 1i*kx_vec .* exp( 1i*kx_vec*dx/2), 1);
ddx_k_shift_neg = ifftshift( 1i*kx_vec .* exp(-1i*kx_vec*dx/2), 1);
ddy_k_shift_pos = ifftshift( 1i*ky_vec .* exp( 1i*ky_vec*dy/2), 2); 
ddy_k_shift_neg = ifftshift( 1i*ky_vec .* exp(-1i*ky_vec*dy/2), 2);
ddz_k_shift_pos = ifftshift( 1i*kz_vec .* exp( 1i*kz_vec*dz/2), 3);
ddz_k_shift_neg = ifftshift( 1i*kz_vec .* exp(-1i*kz_vec*dz/2), 3);
    
% create vector shift operators
x_shift_neg = ifftshift( exp(-1i*kx_vec*dx/2), 1);
y_shift_neg = ifftshift( exp(-1i*ky_vec*dy/2), 2);
z_shift_neg = ifftshift( exp(-1i*kz_vec*dz/2), 3);

% create reduced variables for use with real-to-complex FFT
Nx_r                    = floor(Nx/2) + 1;
Ny_r                    = floor(Ny/2) + 1;
Nz_r                    = floor(Nz/2) + 1;
ddx_k_shift_pos_r       = ddx_k_shift_pos(1:Nx_r);
ddx_k_shift_neg_r       = ddx_k_shift_neg(1:Nx_r);
x_shift_neg_r           = x_shift_neg(1:Nx_r);
y_shift_neg_r           = y_shift_neg(1:Ny_r);
z_shift_neg_r           = z_shift_neg(1:Nz_r);

% create vector PML variables
pml_x       = getPML(Nx, dx, dt, c_ref, pml_x_size, pml_x_alpha, false, 1);
pml_x_sgx   = getPML(Nx, dx, dt, c_ref, pml_x_size, pml_x_alpha, true,  1);
pml_y       = getPML(Ny, dy, dt, c_ref, pml_y_size, pml_y_alpha, false, 2);
pml_y_sgy   = getPML(Ny, dy, dt, c_ref, pml_y_size, pml_y_alpha, true,  2);
pml_z       = getPML(Nz, dz, dt, c_ref, pml_z_size, pml_z_alpha, false, 3);
pml_z_sgz   = getPML(Nz, dz, dt, c_ref, pml_z_size, pml_z_alpha, true,  3);

% cleanup unused variables
clear ddx_k_shift_pos ddx_k_shift_neg x_shift_neg y_shift_neg z_shift_neg;

% =========================================================================
% STORE FLOATS
% =========================================================================

% list of variables stored as floats
variable_names = {...
    'dt', 'dx', 'dy', 'dz', ...
    'ddx_k_shift_pos_r', 'ddx_k_shift_neg_r', 'x_shift_neg_r', ...
    'ddy_k_shift_pos',   'ddy_k_shift_neg',   'y_shift_neg_r', ...
    'ddz_k_shift_pos',   'ddz_k_shift_neg',   'z_shift_neg_r', ...
    'pml_x_sgx',   'pml_y_sgy',   'pml_z_sgz', ...
    'pml_x',       'pml_y',       'pml_z', ...
    'pml_x_alpha', 'pml_y_alpha', 'pml_z_alpha', ...
    'c_ref'};

% change float variables to be in single precision (float in C++), then
% add to HDF5 file
for index = 1:length(variable_names)

    % cast matrix to single precision
    eval([variable_names{index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end

% =========================================================================
% STORE INTEGERS
% =========================================================================

% integer variables
variable_names = {'Nx', 'Ny', 'Nz', 'Nt',...
    'pml_x_size' , 'pml_y_size' , 'pml_z_size'};

% change all the index variables to be in 64-bit unsigned integers (long in C++)
for index = 1:length(variable_names)

    % cast matrix to 64-bit unsigned integer
    eval([variable_names{index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' variable_names{index} ');']);

    % write to HDF5 file
    writeMatrix(filename, eval(variable_names{index}), variable_names{index});

end