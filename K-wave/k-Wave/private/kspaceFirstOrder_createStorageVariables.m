% DESCRIPTION:
%       subscript to create storage variables
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 1st August 2011
%       last update - 8th July 2014
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

% =========================================================================
% PREPARE DATA MASKS AND STORAGE VARIABLES
% =========================================================================

% allocate empty sensor structure (avoids error if only saving the
% final field)
sensor_data = struct;

% check sensor mask based on the Cartesian interpolation setting
if ~record.binary_sensor_mask && strcmp(cartesian_interp, 'nearest')

    % extract the data using the binary sensor mask created in
    % inputChecking, but switch on Cartesian reorder flag so that the
    % final data is returned in the correct order (not in time
    % reversal mode).
    record.binary_sensor_mask = true;
    if ~record.time_rev
        record.reorder_data = true;
    end

    % check if any duplicate points have been discarded in the
    % conversion from a Cartesian to binary mask
    num_discarded_points = length(sensor_x) - sum(sensor.mask(:));
    if num_discarded_points ~= 0
        disp(['  WARNING: ' num2str(num_discarded_points) ' duplicated sensor points discarded (nearest neighbour interpolation)']);
    end        

end

% preallocate output variables
if ~record.time_rev

    % get the number of sensor points
    if record.blank_sensor
        num_sensor_points = numel(kgrid.k);
    elseif record.binary_sensor_mask
        num_sensor_points = length(sensor_mask_index);
    else
        num_sensor_points = length(sensor_x);
    end

    % calculate the number of time points that are stored.
    % - if streaming data to disk, reduce to the size of the
    %   sensor_data matrix based on the value of record.stream_to_disk
    % - if a user input for sensor.record_start_index is given, reduce
    %   the size of the sensor_data matrix based on the value given
    if kgrid.dim == 3 && record.stream_to_disk
        num_recorded_time_points = record.stream_to_disk;

        % initialise the file index variable
        stream_data_index = 1;
    else
        num_recorded_time_points = kgrid.Nt - sensor.record_start_index + 1;
    end

    % create shift operators used for calculating the components of the
    % particle velocity field on the non-staggered grids (these are
    % used for both binary and cartesian sensor masks)
    if (record.u_non_staggered || record.I || record.I_avg)
        if use_sg
            switch kgrid.dim
                case 1
                    record.x_shift_neg = ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
                case 2
                    record.x_shift_neg = ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
                    record.y_shift_neg = ifftshift( exp(-1i*kgrid.ky_vec*kgrid.dy/2) ).';
                case 3
                    record.x_shift_neg = ifftshift( exp(-1i*kgrid.kx_vec*kgrid.dx/2) );
                    record.y_shift_neg = ifftshift( exp(-1i*kgrid.ky_vec*kgrid.dy/2) ).';
                    record.z_shift_neg = permute(ifftshift( exp(-1i*kgrid.kz_vec*kgrid.dz/2) ), [2 3 1]);
            end
        else
            switch kgrid.dim
                case 1
                    record.x_shift_neg = 1;
                case 2
                    record.x_shift_neg = 1;
                    record.y_shift_neg = 1;
                case 3
                    record.x_shift_neg = 1;
                    record.y_shift_neg = 1;
                    record.z_shift_neg = 1;
            end
        end
    end

    % create storage and scaling variables - all variables are saved as
    % fields of a structure called sensor_data. If only p is being stored
    % (i.e., if no user input is given for sensor.record), then
    % sensor_data.p is copied to sensor_data at the end of the simulation

    % time history of the acoustic pressure
    if record.p || record.I || record.I_avg
        sensor_data.p = castZeros([num_sensor_points, num_recorded_time_points]);
    end

    % maximum pressure
    if record.p_max
        sensor_data.p_max = castZeros([num_sensor_points, 1]);
    end

    % minimum pressure
    if record.p_min
        sensor_data.p_min = castZeros([num_sensor_points, 1]);
    end        

    % rms pressure
    if record.p_rms
        sensor_data.p_rms = castZeros([num_sensor_points, 1]);
    end

    % calculate the size of the _all and _final output variables. If
    % the PML is set to be outside the grid, these will be the same
    % size as the user input, rather than the expanded grid
    if PML_inside
        all_vars_size = size(kgrid.k);
    else
        switch kgrid.dim
            case 1
                all_vars_size = [kgrid.Nx - 2*PML_x_size, 1];
            case 2
                all_vars_size = [kgrid.Nx - 2*PML_x_size, kgrid.Ny - 2*PML_y_size];
            case 3
                all_vars_size = [kgrid.Nx - 2*PML_x_size, kgrid.Ny - 2*PML_y_size, kgrid.Nz - 2*PML_z_size];
        end
    end

    % maximum pressure over all grid points
    if record.p_max_all
        sensor_data.p_max_all = castZeros(all_vars_size);
    end

    % minimum pressure over all grid points
    if record.p_min_all
        sensor_data.p_min_all = castZeros(all_vars_size);
    end           

    % time history of the acoustic particle velocity
    if record.u
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1
                sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
            case 2
                sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uy = castZeros([num_sensor_points, num_recorded_time_points]);
            case 3
                sensor_data.ux = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uy = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uz = castZeros([num_sensor_points, num_recorded_time_points]); 
        end                
    end        

    % maximum particle velocity
    if record.u_max
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1            
                sensor_data.ux_max = castZeros([num_sensor_points, 1]);
            case 2
                sensor_data.ux_max = castZeros([num_sensor_points, 1]);
                sensor_data.uy_max = castZeros([num_sensor_points, 1]);
            case 3
                sensor_data.ux_max = castZeros([num_sensor_points, 1]);
                sensor_data.uy_max = castZeros([num_sensor_points, 1]);
                sensor_data.uz_max = castZeros([num_sensor_points, 1]);
        end
    end

    % minimum particle velocity
    if record.u_min
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1            
                sensor_data.ux_min = castZeros([num_sensor_points, 1]);
            case 2
                sensor_data.ux_min = castZeros([num_sensor_points, 1]);
                sensor_data.uy_min = castZeros([num_sensor_points, 1]);
            case 3
                sensor_data.ux_min = castZeros([num_sensor_points, 1]);
                sensor_data.uy_min = castZeros([num_sensor_points, 1]);
                sensor_data.uz_min = castZeros([num_sensor_points, 1]);
        end
    end        

    % rms particle velocity
    if record.u_rms
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1            
                sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
            case 2
                sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
                sensor_data.uy_rms = castZeros([num_sensor_points, 1]);
            case 3
                sensor_data.ux_rms = castZeros([num_sensor_points, 1]);
                sensor_data.uy_rms = castZeros([num_sensor_points, 1]);
                sensor_data.uz_rms = castZeros([num_sensor_points, 1]);
        end
    end

    % maximum particle velocity over all grid points
    if record.u_max_all
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1            
                sensor_data.ux_max_all = castZeros(all_vars_size);
            case 2
                sensor_data.ux_max_all = castZeros(all_vars_size);
                sensor_data.uy_max_all = castZeros(all_vars_size);
            case 3
                sensor_data.ux_max_all = castZeros(all_vars_size);
                sensor_data.uy_max_all = castZeros(all_vars_size);
                sensor_data.uz_max_all = castZeros(all_vars_size);
        end
    end

    % minimum particle velocity over all grid points
    if record.u_min_all
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1            
                sensor_data.ux_min_all = castZeros(all_vars_size);
            case 2
                sensor_data.ux_min_all = castZeros(all_vars_size);
                sensor_data.uy_min_all = castZeros(all_vars_size);
            case 3
                sensor_data.ux_min_all = castZeros(all_vars_size);
                sensor_data.uy_min_all = castZeros(all_vars_size);
                sensor_data.uz_min_all = castZeros(all_vars_size);
        end
    end           

    % time history of the acoustic particle velocity on the
    % non-staggered grid points
    if record.u_non_staggered || record.I || record.I_avg
        % pre-allocate the velocity fields based on the number of
        % dimensions in the simulation
        switch kgrid.dim
            case 1
                sensor_data.ux_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]);
            case 2
                sensor_data.ux_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uy_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]);
            case 3
                sensor_data.ux_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uy_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]);
                sensor_data.uz_non_staggered = castZeros([num_sensor_points, num_recorded_time_points]); 
        end                
    end        

    % object of the kWaveTransducer class is being used as a sensor            
    if record.transducer_sensor
        if record.transducer_receive_elevation_focus
            % if there is elevation focusing, a buffer is
            % needed to store a short time history at each
            % sensor point before averaging
            sensor_data_buffer_size = max(sensor.elevation_beamforming_delays) + 1;
            if sensor_data_buffer_size > 1
                sensor_data_buffer = castZeros([num_sensor_points, sensor_data_buffer_size]);
            else
                clear sensor_data_buffer sensor_data_buffer_size;
                record.transducer_receive_elevation_focus = false;
            end
        end

        % the grid points can be summed on the fly and so the
        % sensor is the size of the number of active elements 
        sensor_data.transducer = castZeros([sensor.number_active_elements, num_recorded_time_points]);
    end

    % precomputate the triangulation points if a Cartesian sensor mask
    % is used with linear interpolation (tri and bc are the Delaunay
    % triangulation and Barycentric coordinates)
    if ~record.binary_sensor_mask
        if kgrid.dim == 1

            % assign pseudonym for Cartesain grid points in 1D (this
            % is later used for data casting)
            record.grid_x = kgrid.x_vec;

        else

            % update command line status
            disp('  calculating Delaunay triangulation...');

            % compute triangulation
            if kgrid.dim == 2
                [record.tri, record.bc] = gridDataFast2D(kgrid.x, kgrid.y, sensor_x, sensor_y);
            elseif kgrid.dim == 3
                [record.tri, record.bc] = gridDataFast3D(kgrid.x, kgrid.y, kgrid.z, sensor_x, sensor_y, sensor_z);
            end

        end
    end
end