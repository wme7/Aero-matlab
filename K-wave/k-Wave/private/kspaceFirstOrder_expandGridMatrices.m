% DESCRIPTION:
%       script to expand the grid matrices used in kspaceFirstOrder1D,
%       kspaceFirstOrder2D, and kspaceFirstOrder3D
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 20th August 2010
%       last update - 25th August 2014
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

%#ok<*STISA>

% update command line status
disp('  expanding computational grid...');

% retaining the values for kgrid.t_array
t_array_temp = kgrid.t_array;

% expand the computational grid, replacing the original grid
switch kgrid.dim
    case 1
        kgrid = makeGrid(kgrid.Nx + 2*PML_x_size, kgrid.dx);
    case 2
        kgrid = makeGrid(kgrid.Nx + 2*PML_x_size, kgrid.dx, kgrid.Ny + 2*PML_y_size, kgrid.dy);
    case 3
        kgrid = makeGrid(kgrid.Nx + 2*PML_x_size, kgrid.dx, kgrid.Ny + 2*PML_y_size, kgrid.dy, kgrid.Nz + 2*PML_z_size, kgrid.dz);
end

% re-assign original time array
kgrid.t_array = t_array_temp;
clear t_array_temp;

% assign Nt and dt to kgrid if given as a structure
if isstruct(kgrid)
    kgrid.Nt = length(kgrid.t_array);
    kgrid.dt = kgrid.t_array(2) - kgrid.t_array(1);
end 

% set the PML size for use with expandMatrix
switch kgrid.dim
    case 1
        expand_size = PML_x_size;
    case 2
        expand_size = [PML_x_size, PML_y_size]; 
    case 3
        expand_size = [PML_x_size, PML_y_size, PML_z_size];
end

% update the data type in case adding the PML requires additional index
% precision
if kgrid.total_grid_points < intmax('uint8');
    index_data_type = 'uint8';
elseif kgrid.total_grid_points < intmax('uint16');
    index_data_type = 'uint16';
elseif kgrid.total_grid_points < intmax('uint32');
    index_data_type = 'uint32';                
else
    index_data_type = 'double';
end  

% enlarge the sensor mask (for Cartesian sensor masks and cuboid corners,
% this has already been converted to a binary mask for display in
% inputChecking)
if record.use_sensor && ~record.blank_sensor
    if strcmp(class(sensor), 'kWaveTransducer') 
        sensor.expand_grid(expand_size);
    else
        sensor.mask = expandMatrix(sensor.mask, expand_size, 0);
    end 
end

% add the PML size to cuboid corner indices if using a cuboid sensor mask
if record.cuboid_corners
   switch kgrid.dim
        case 1
            record.cuboid_corners_list = record.cuboid_corners_list + PML_x_size;
        case 2
            record.cuboid_corners_list([1, 3], :) = record.cuboid_corners_list([1, 3], :) + PML_x_size;
            record.cuboid_corners_list([2, 4], :) = record.cuboid_corners_list([2, 4], :) + PML_y_size;
        case 3
            record.cuboid_corners_list([1, 4], :) = record.cuboid_corners_list([1, 4], :) + PML_x_size;
            record.cuboid_corners_list([2, 5], :) = record.cuboid_corners_list([2, 5], :) + PML_y_size;            
            record.cuboid_corners_list([3, 6], :) = record.cuboid_corners_list([3, 6], :) + PML_z_size;
    end 
end

% enlarge the sound speed grids by extending the edge values into the
% expanded grid
if ~elastic_code
    if numel(medium.sound_speed) > 1
        medium.sound_speed = expandMatrix(medium.sound_speed, expand_size);
    end   
else
    if numel(medium.sound_speed_compression) > 1
        medium.sound_speed_compression = expandMatrix(medium.sound_speed_compression, expand_size);
    end
    if numel(medium.sound_speed_shear) > 1
        medium.sound_speed_shear = expandMatrix(medium.sound_speed_shear, expand_size);
    end
end

% enlarge the grid of density by extending the edge values into the
% expanded grid    
if numel(medium.density) > 1
    medium.density = expandMatrix(medium.density, expand_size);
end

% enlarge the grid of medium.alpha_coeff if given
if isfield(medium, 'alpha_coeff') && numel(medium.alpha_coeff) > 1
    medium.alpha_coeff = expandMatrix(medium.alpha_coeff, expand_size);
end

% enlarge the grid of medium.alpha_coeff_compression if given
if isfield(medium, 'alpha_coeff_compression') && numel(medium.alpha_coeff_compression) > 1
    medium.alpha_coeff_compression = expandMatrix(medium.alpha_coeff_compression, expand_size);
end

% enlarge the grid of medium.alpha_coeff_shear if given
if isfield(medium, 'alpha_coeff_shear') && numel(medium.alpha_coeff_shear) > 1
    medium.alpha_coeff_shear = expandMatrix(medium.alpha_coeff_shear, expand_size);
end

% enlarge the grid of medium.BonA if given
if isfield(medium, 'BonA') && numel(medium.BonA) > 1
    medium.BonA = expandMatrix(medium.BonA, expand_size);
end

% enlarge the display mask if given
if ~(strcmp(display_mask, 'default') || strcmp(display_mask, 'off'))
    display_mask = expandMatrix(display_mask, expand_size, 0);
end

% enlarge the initial pressure if given
if isfield(source, 'p0')
    source.p0 = expandMatrix(source.p0, expand_size, 0);
end      

% enlarge the absorption filter mask if given
if isfield(medium, 'alpha_filter');
    medium.alpha_filter = expandMatrix(medium.alpha_filter, expand_size, 0);
end 

% enlarge the pressure source mask if given
if p_source   
    
    % enlarge the pressure source mask
    source.p_mask = expandMatrix(source.p_mask, expand_size, 0);

    % create an indexing variable corresponding to the source elements
    p_source_index = find(source.p_mask ~= 0);
    
    % convert the data type depending on the number of indices
    eval(['p_source_index = ' index_data_type '(p_source_index);']);      

end

% enlarge the velocity source mask if given
if (ux_source || uy_source || uz_source || transducer_source)
    
    % update the source indexing variable
    if strcmp(class(source), 'kWaveTransducer')
        
        % check if the sensor is also the same transducer, if so, don't
        % expand the grid again
        if ~(strcmp(class(sensor), 'kWaveTransducer') && isequal(sensor, source))
        
            % expand the transducer mask
            source.expand_grid(expand_size);
            
        end
        
        % get the new active elements mask
        active_elements_mask = source.active_elements_mask;

        % update the indexing variable corresponding to the active elements
        u_source_index = find(active_elements_mask ~= 0);
        
        % clean up unused variables
        clear active_elements_mask;        
        
    else
        
        % enlarge the velocity source mask
        source.u_mask = expandMatrix(source.u_mask, expand_size, 0);
        
        % create an indexing variable corresponding to the source elements
        u_source_index = find(source.u_mask ~= 0);  
        
    end
    
    % convert the data type depending on the number of indices
    eval(['u_source_index = ' index_data_type '(u_source_index);']);   
    
end

% enlarge the stress source mask if given
if (sxx_source || syy_source || szz_source || sxy_source || sxz_source || syz_source)
   
    % enlarge the velocity source mask
    source.s_mask = expandMatrix(source.s_mask, expand_size, 0);

    % create an indexing variable corresponding to the source elements
    s_source_index = find(source.s_mask ~= 0);  
    
    % convert the data type depending on the number of indices
    eval(['s_source_index = ' index_data_type '(s_source_index);']);  
    
end

% enlarge the directivity angle if given (2D only)
if record.use_sensor && kgrid.dim == 2 && record.compute_directivity
    
    % enlarge the directivity angle
    sensor.directivity_angle = expandMatrix(sensor.directivity_angle, expand_size, 0);
    
    % re-assign the wavenumber vectors
    sensor.directivity_wavenumbers = [kgrid.ky(:)'; kgrid.kx(:)'];
    
end    

% update command line status
switch kgrid.dim
    case 1
        disp(['  computation grid size: ' num2str(kgrid.Nx) ' grid points']);
    case 2
        disp(['  computational grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' grid points']);
    case 3
        disp(['  computational grid size: ' num2str(kgrid.Nx) ' by ' num2str(kgrid.Ny) ' by ' num2str(kgrid.Nz) ' grid points']);
end