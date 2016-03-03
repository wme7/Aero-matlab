% DESCRIPTION:
%       subscript to save input data to disk
%
% ABOUT:
%       author      - Bradley Treeby and Jiri Jaros
%       date        - 24th August 2011
%       last update - 21st August 2014
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

% update command line status
disp(['  precomputation completed in ' scaleTime(toc)]);
tic;
disp('  saving input files to disk...');

% =========================================================================
% VARIABLE LIST - THESE ARE USED IN ALL SIMULATIONS
% =========================================================================

% list of all the single precision variables used within the time loop in
% both the fluid and elastic codes
variable_list = {'dt', 'c_ref', ...
    'ddx_k_shift_pos_r', 'ddx_k_shift_neg_r', 'x_shift_neg_r', ...
    'ddy_k_shift_pos',   'ddy_k_shift_neg',   'y_shift_neg_r', ... 
    'ddz_k_shift_pos',   'ddz_k_shift_neg',   'z_shift_neg_r', ...
    };

% list of all the integer variables used within the time loop regardless of
% options
integer_variable_list = {...
    'Nx', 'Ny', 'Nz', 'Nt', ...
    'ux_source_flag', ...
    'uy_source_flag', ...
    'uz_source_flag', ...
    'p_source_flag', ...
    'p0_source_flag', ...
    'transducer_source_flag', ...
    'nonuniform_grid_flag', ...
    'nonlinear_flag', ...
    'absorbing_flag', ...
    'sensor_mask_type', ...
    };

% additional single precision variables not used within time loop but
% stored directly to output file
variable_list = [variable_list, {...
    'dx', 'dy', 'dz'...
    'pml_x_alpha', 'pml_y_alpha', 'pml_z_alpha'...
    }];

% additional integer variables not used within time loop but stored
% directly to output file
integer_variable_list = [integer_variable_list, {...
    'pml_x_size', 'pml_y_size', 'pml_z_size'...
    }];

% =========================================================================
% GRID VARIABLES
% =========================================================================

% create pseudonyms for variables stored in structures
Nx                      = kgrid.Nx;
Ny                      = kgrid.Ny;
Nz                      = kgrid.Nz;
Nt                      = length(t_array);
dx                      = kgrid.dx;
dy                      = kgrid.dy;
dz                      = kgrid.dz;

% create shift variables used for calculating u_non_staggered and I outputs
x_shift_neg             =         ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
y_shift_neg             =         ifftshift( exp(-1i*kgrid.ky_vec*kgrid.dy/2) ).';
z_shift_neg             = permute(ifftshift( exp(-1i*kgrid.kz_vec*kgrid.dz/2) ), [2 3 1]);

% create reduced variables for use with real-to-complex FFT
Nx_r                    = floor(Nx/2) + 1;
Ny_r                    = floor(Ny/2) + 1;
Nz_r                    = floor(Nz/2) + 1;
ddx_k_shift_pos_r       = ddx_k_shift_pos(1:Nx_r);
ddx_k_shift_neg_r       = ddx_k_shift_neg(1:Nx_r);
x_shift_neg_r           = x_shift_neg(1:Nx_r);
y_shift_neg_r           = y_shift_neg(1:Ny_r);
z_shift_neg_r           = z_shift_neg(1:Nz_r);

% create pseudonyms for variables that have different names to the MATLAB
% version
pml_x_size              = PML_x_size;
pml_y_size              = PML_y_size;
pml_z_size              = PML_z_size;
pml_x_alpha             = PML_x_alpha;
pml_y_alpha             = PML_y_alpha;
pml_z_alpha             = PML_z_alpha;
c0                      = c;

% create pseudonyms for the source flags
ux_source_flag          = ux_source;
uy_source_flag          = uy_source;
uz_source_flag          = uz_source;
p_source_flag           = p_source;
p0_source_flag          = isfield(source, 'p0');
elastic_flag            = elastic_code;
transducer_source_flag  = transducer_source;

% create pseudonyms for the sensor flags
%   0: binary mask indices
%   1: cuboid corners
sensor_mask_type        = record.cuboid_corners;

% cleanup unused variables
clear x_shift_neg;

% =========================================================================
% MEDIUM PROPERTIES
% =========================================================================

if elastic_code
    
    % add elastic medium properties to the variable list
    variable_list = [variable_list, {'lambda',...
        'mu', 'mu_sgxy', 'mu_sgxz', 'mu_sgyz', ...
        'rho0_sgx', 'rho0_sgy', 'rho0_sgz'}];
    
else
    
    % add fluid medium properties to the variable list
    variable_list = [variable_list, {'c0', ...
        'rho0', 'rho0_sgx', 'rho0_sgy', 'rho0_sgz'}];
    
end

% =========================================================================
% PML PROPERTIES
% =========================================================================

if elastic_code
    
    % add elastic PML variables to the list
    variable_list = [variable_list, {...
        'pml_x_sgx',  'pml_y_sgy',  'pml_z_sgz',...
        'pml_x',      'pml_y',      'pml_z',...
        'mpml_x_sgx', 'mpml_y_sgy', 'mpml_z_sgz',...
        'mpml_x',     'mpml_y',     'mpml_z'}];
    
else
    
    % add fluid PML variables to the list
    variable_list = [variable_list, {...
        'pml_x_sgx', 'pml_y_sgy', 'pml_z_sgz',...
        'pml_x',     'pml_y',     'pml_z'}];
    
end

% =========================================================================
% VARIABLES USED IN NONLINEAR SIMULATIONS
% =========================================================================

if nonlinear
    
    % set nonlinear flag
    nonlinear_flag = 1;
    
    % create pseudonyms for variables saved in structures
    BonA = medium.BonA;
    
    % add BonA to the variable list
    variable_list = [variable_list, {'BonA'}];

else
    
    % set nonlinear flag
    nonlinear_flag = 0;
    
end

% =========================================================================
% VARIABLES USED IN ABSORBING SIMULATIONS
% =========================================================================

if strcmp(equation_of_state, 'absorbing')
    
    % set absorbing flag
    absorbing_flag = 1;  
    
    % create pseudonyms for absorption variables saved in structures
    alpha_coeff = medium.alpha_coeff;
    alpha_power = medium.alpha_power;
    
    % add to the variable list
    variable_list = [variable_list, {'alpha_coeff', 'alpha_power'}];

else
    
    % set absorbing flag
    absorbing_flag = 0;  
    
end
    
% =========================================================================
% SOURCE VARIABLES
% =========================================================================
% source modes and indicies
% - these are only defined if the source flags are > 0
% - the source mode describes whether the source will be added or replaced
% - the source indicies describe which grid points act as the source
% - the u_source_index is reused for any of the u sources and the transducer source

% velocity sources
if ux_source_flag || uy_source_flag || uz_source_flag
    u_source_mode = ~strcmp(source.u_mode, 'dirichlet');
    if ux_source_flag
        u_source_many = numDim(source.ux) > 1;
    elseif uy_source_flag
        u_source_many = numDim(source.uy) > 1;
    elseif uz_source_flag
        u_source_many = numDim(source.uz) > 1;
    end
    integer_variable_list = [integer_variable_list, {'u_source_mode', 'u_source_many', 'u_source_index'}];
end

% pressure source
if p_source_flag
    p_source_mode = ~strcmp(source.p_mode, 'dirichlet');
    p_source_many = numDim(source.p) > 1;
    integer_variable_list = [integer_variable_list, {'p_source_mode', 'p_source_many', 'p_source_index'}];
end

% transducer source
if transducer_source_flag
    integer_variable_list = [integer_variable_list, {'u_source_index'}];
end

% source variables
% - these are only defined if the source flags are > 0
% - these are the actual source values
% - these are indexed as (position_index, time_index)
if ux_source_flag
    ux_source_input = source.ux;
    variable_list = [variable_list, {'ux_source_input'}];
end
if uy_source_flag
    uy_source_input = source.uy;
    variable_list = [variable_list, {'uy_source_input'}];
end
if uz_source_flag 
    uz_source_input = source.uz;
    variable_list = [variable_list, {'uz_source_input'}];
end
if p_source_flag
    p_source_input = source.p;
    variable_list = [variable_list, {'p_source_input'}];
end
if transducer_source_flag
    transducer_source_input = transducer_input_signal;
    variable_list = [variable_list, {'transducer_source_input'}];
    integer_variable_list = [integer_variable_list, {'delay_mask'}];
end

% initial pressure source variable
% - this is only defined if the p0 source flag is 1
% - this defines the initial pressure everywhere (there is no indicies)
if p0_source_flag
    p0_source_input = source.p0;
    variable_list = [variable_list, {'p0_source_input'}];
end

% =========================================================================
% SENSOR VARIABLES
% =========================================================================

if sensor_mask_type == 0
    % mask is defined as a list of grid indices
    integer_variable_list = [integer_variable_list, {'sensor_mask_index'}];
elseif sensor_mask_type == 1
    % mask is defined as a list of cuboid corners
    sensor_mask_corners = record.cuboid_corners_list;
    integer_variable_list = [integer_variable_list, {'sensor_mask_corners'}];
else
    error('unknown option for sensor_mask_type');
end

% =========================================================================
% VARIABLES USED FOR NONUNIFORM GRIDS
% =========================================================================

% set nonuniform flag and variables
% - these are only defined if nonuniform_grid_flag is 1
% - these are applied using the bsxfun formulation
nonuniform_grid_flag = nonuniform_grid;
if nonuniform_grid_flag
    variable_list = [variable_list, {'dxudxn', 'dyudyn', 'dzudzn', 'dxudxn_sgx', 'dyudyn_sgy', 'dzudzn_sgz'}];
    dxudxn = kgrid.dxudxn;
    if numel(dxudxn) == 1
        dxudxn = ones(kgrid.Nx, 1);
    end
    dyudyn = kgrid.dyudyn;
    if numel(dyudyn) == 1
        dyudyn = ones(1, kgrid.Ny);
    end
    dzudzn = kgrid.dzudzn;
    if numel(dzudzn) == 1
        dzudzn = ones(1, 1, kgrid.Nz);
    end
    dxudxn_sgx = kgrid.dxudxn_sgx;
    if numel(dxudxn) == 1
        dxudxn_sgx = ones(kgrid.Nx, 1);
    end
    dyudyn_sgy = kgrid.dyudyn_sgy;
    if numel(dyudyn) == 1
        dyudyn_sgy = ones(1, kgrid.Ny);
    end
    dzudzn_sgz = kgrid.dzudzn_sgz;
    if numel(dzudzn) == 1
        dzudzn_sgz = ones(1, 1, kgrid.Nz);
    end    
end

% =========================================================================
% DATACAST AND SAVING
% =========================================================================

% check for HDF5 filename extension
[~, ~, filename_ext] = fileparts(save_to_disk);

% use .h5 as default if no extension is given
if isempty(filename_ext)
    filename_ext = '.h5';
    save_to_disk = [save_to_disk '.h5'];
end

% save file
if strcmp(filename_ext, '.h5')
    
    % ----------------
    % SAVE HDF5 FILE
    % ----------------
    
    % check if file exists, and delete if it does (the hdf5 library will
    % give an error if the file already exists)
    if exist(save_to_disk, 'file')
        delete(save_to_disk);
    end

    % get HDF5 literals
    getH5Literals;

    % change all the variables to be in single precision (float in C++),
    % then add to HDF5 File
    for cast_index = 1:length(variable_list)

        % cast matrix to single precision
        eval([variable_list{cast_index} ' = ' MATRIX_DATA_TYPE_MATLAB '(' variable_list{cast_index} ');']);

        % write to HDF5 file
        writeMatrix(save_to_disk, eval(variable_list{cast_index}), variable_list{cast_index});
        
        % clear matrix from memory
        eval(['clear ' variable_list{cast_index}]);

    end

    % change all the index variables to be in 64-bit unsigned integers
    % (long in C++), then add to HDF5 file
    for cast_index = 1:length(integer_variable_list)

        % cast matrix to 64-bit unsigned integer
        eval([integer_variable_list{cast_index} ' = ' INTEGER_DATA_TYPE_MATLAB '(' integer_variable_list{cast_index} ');']);

        % write to HDF5 file
        writeMatrix(save_to_disk, eval(integer_variable_list{cast_index}), integer_variable_list{cast_index});

        % clear matrix from memory
        eval(['clear ' integer_variable_list{cast_index}]);
        
    end

    % set additional file attributes
    writeAttributes(save_to_disk);
    
elseif strcmp(filename_ext, '.mat')

    % ----------------
    % SAVE .MAT FILE
    % ----------------

    % change all the variables to be in single precision (float in C++)
    data_cast = 'single';
    for cast_index = 1:length(variable_list)
        eval([variable_list{cast_index} ' = ' data_cast '(' variable_list{cast_index} ');']);
    end

    % change all the index variables to be in 64-bit unsigned integers (long in C++)
    data_cast = 'uint64';
    for cast_index = 1:length(integer_variable_list)
        eval([integer_variable_list{cast_index} ' = ' data_cast '(' integer_variable_list{cast_index} ');']);
    end

    % save the input variables to disk as a MATLAB binary file
    variable_list = [variable_list, integer_variable_list];
    save(save_to_disk, variable_list{:}, '-v7.3');   
    
else
    
    % throw error for unknown filetype
    error('unknown file extension for ''SaveToDisk'' filename');
    
end

% update command line status
disp(['  completed in ' scaleTime(toc)]);