% DESCRIPTION:
%       Subscript to calculate the acoustic intensity from the time varying
%       acoustic pressure and particle velocity recorded during the
%       simulation. The particle velocity is first temporally shifted
%       forwards by dt/2 using a spectral method so both variables are on
%       the regular (non-staggered) grid.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 4th September 2013
%       last update - 28th August 2014
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

%#ok<*SAGROW>

% check if using cuboid corners or a binary/cartesian sensor mask
if save_intensity_matlab_code || ~record.cuboid_corners
    
    % --------------------------------------------------------------------
    % SENSOR MASK FROM MATLAB CODE ORDERED AS (SENSOR_INDEX, TIME_INDEX)
    % --------------------------------------------------------------------
    
    % shift the recorded particle velocity to the regular
    % (non-staggered) temporal grid
    ux_sgt = timeShift(sensor_data.ux_non_staggered, kgrid.dt);
    if kgrid.dim > 1
        uy_sgt = timeShift(sensor_data.uy_non_staggered, kgrid.dt);
    end
    if kgrid.dim > 2        
        uz_sgt = timeShift(sensor_data.uz_non_staggered, kgrid.dt);
    end 
    
    % compute the time varying intensity
    sensor_data.Ix = sensor_data.p .* ux_sgt;
    if kgrid.dim > 1
        sensor_data.Iy = sensor_data.p .* uy_sgt;
    end
    if kgrid.dim > 2
        sensor_data.Iz = sensor_data.p .* uz_sgt;
    end  

    % calculate the time average of the intensity if required using the last
    % dimension (this works for both linear and cuboid sensor masks)
    if record.I_avg
        sensor_data.Ix_avg = mean(sensor_data.Ix, 2);
        if kgrid.dim > 1
            sensor_data.Iy_avg = mean(sensor_data.Iy, 2);
        end
        if kgrid.dim > 2
            sensor_data.Iz_avg = mean(sensor_data.Iz, 2);
        end
    end
    
else
    
    % --------------------------------------------------------------------
    % CUBOID CORNERS SENSOR MASK FROM C++ CODE 
    % --------------------------------------------------------------------
    
    % loop over the number of cuboid corners
    for cuboid_index = 1:size(record.cuboid_corners_list, 2)
    
        % extract the size of the time varying outputs for the current
        % cuboid
        output_size = size(sensor_data(cuboid_index).ux_non_staggered);
        
        % make a copy of the staggered grid velocity for shifting
        ux_sgt = sensor_data(cuboid_index).ux_non_staggered;
        if kgrid.dim > 1
            uy_sgt = sensor_data(cuboid_index).uy_non_staggered;
        end
        if kgrid.dim > 2        
            uz_sgt = sensor_data(cuboid_index).uz_non_staggered;
        end 

        % reshape to (index, time)
        ux_sgt = reshape(ux_sgt, [], output_size(end)); 
        if kgrid.dim > 1
            uy_sgt = reshape(uy_sgt, [], output_size(end));
        end
        if kgrid.dim > 2        
            uz_sgt = reshape(uz_sgt, [], output_size(end));
        end         
        
        % shift the recorded particle velocity to the regular
        % (non-staggered) temporal grid
        ux_sgt = timeShift(ux_sgt, kgrid.dt);
        if kgrid.dim > 1
            uy_sgt = timeShift(uy_sgt, kgrid.dt);
        end
        if kgrid.dim > 2        
            uz_sgt = timeShift(uz_sgt, kgrid.dt);
        end 

        % reshape to original size
        ux_sgt = reshape(ux_sgt, output_size);
        if kgrid.dim > 1
            uy_sgt = reshape(uy_sgt, output_size);
        end
        if kgrid.dim > 2        
            uz_sgt = reshape(uz_sgt, output_size);
        end 
        
        % compute the time varying intensity
        sensor_data(cuboid_index).Ix = sensor_data(cuboid_index).p .* ux_sgt;
        if kgrid.dim > 1
            sensor_data(cuboid_index).Iy = sensor_data(cuboid_index).p .* uy_sgt;
        end
        if kgrid.dim > 2
            sensor_data(cuboid_index).Iz = sensor_data(cuboid_index).p .* uz_sgt;
        end  

        % calculate the time average of the intensity if required using the last
        % dimension
        if record.I_avg
            sensor_data(cuboid_index).Ix_avg = mean(sensor_data(cuboid_index).Ix, numel(output_size));
            if kgrid.dim > 1
                sensor_data(cuboid_index).Iy_avg = mean(sensor_data(cuboid_index).Iy, numel(output_size));
            end
            if kgrid.dim > 2
                sensor_data(cuboid_index).Iz_avg = mean(sensor_data(cuboid_index).Iz, numel(output_size));
            end
        end
        
    end
    
end

% remove the non staggered particle velocity variables if not required
if ~record.u_non_staggered
    switch kgrid.dim
        case 1
            sensor_data = rmfield(sensor_data, {'ux_non_staggered'});
        case 2
            sensor_data = rmfield(sensor_data, {'ux_non_staggered', 'uy_non_staggered'});
        case 3
            sensor_data = rmfield(sensor_data, {'ux_non_staggered', 'uy_non_staggered', 'uz_non_staggered'});
    end
end

% remove the time varying intensity if not required
if ~record.I
    switch kgrid.dim
        case 1
            sensor_data = rmfield(sensor_data, {'Ix'});
        case 2
            sensor_data = rmfield(sensor_data, {'Ix', 'Iy'});
        case 3
            sensor_data = rmfield(sensor_data, {'Ix', 'Iy', 'Iz'});
    end
end

% remove the time varying pressure if not required
if ~record.p
    sensor_data = rmfield(sensor_data, {'p'});
end

% remove unused variables
clear save_intensity_matlab_code output_size ux_sgt uy_sgt uz_sgt;

%     % plot an example of the shifted velocity
%     figure;
%     plot((1:length(sensor_data.ux))*kgrid.dt, sensor_data.ux(1, :), 'k-s'); 
%     hold on;
%     plot((1:length(sensor_data.ux))*kgrid.dt + kgrid.dt/2, ux_sgt(1, :), 'r-s');
%     legend('Original', 'Shifted');