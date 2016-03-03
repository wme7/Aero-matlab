% DESCRIPTION:
%       subscript reorder the sensor points if a binary sensor mask was
%       used for Cartesian sensor mask nearest neighbour interpolation
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 2nd September 2012
%       last update - 18th August 2014
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
disp('  reordering Cartesian measurement data...');

if record.p
    % calculate the position for reordering data
    new_col_pos = length(sensor_data.p(1,:)) + 1;
    
    % append the reordering data
    sensor_data.p(:, new_col_pos) = reorder_index;

    % reorder based on the order_index
    sensor_data.p = sortrows(sensor_data.p, new_col_pos);

    % remove the reordering data
    sensor_data.p = sensor_data.p(:, 1:new_col_pos - 1);
end

if record.p_max
    % append the reordering data
    sensor_data.p_max(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.p_max = sortrows(sensor_data.p_max, 2);

    % remove the reordering data
    sensor_data.p_max = sensor_data.p_max(:, 1);
end

if record.p_min
    % append the reordering data
    sensor_data.p_min(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.p_min = sortrows(sensor_data.p_min, 2);

    % remove the reordering data
    sensor_data.p_min = sensor_data.p_min(:, 1);
end

if record.p_rms
    % append the reordering data
    sensor_data.p_rms(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.p_rms = sortrows(sensor_data.p_rms, 2);

    % remove the reordering data
    sensor_data.p_rms = sensor_data.p_rms(:, 1);
end

if record.u
    % calculate the position for reordering data
    new_col_pos = length(sensor_data.ux(1,:)) + 1;    
    
    % append the reordering data
    sensor_data.ux(:, new_col_pos) = reorder_index;

    % reorder based on the order_index
    sensor_data.ux = sortrows(sensor_data.ux, new_col_pos);
    
    % remove the reordering data
    sensor_data.ux = sensor_data.ux(:, 1:new_col_pos - 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy(:, new_col_pos) = reorder_index;       
        sensor_data.uy = sortrows(sensor_data.uy, new_col_pos);
        sensor_data.uy = sensor_data.uy(:, 1:new_col_pos - 1);
    end
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz(:, new_col_pos) = reorder_index; 
        sensor_data.uz = sortrows(sensor_data.uz, new_col_pos);
        sensor_data.uz = sensor_data.uz(:, 1:new_col_pos - 1);           
    end
end

if record.u_max
    % append the reordering data
    sensor_data.ux_max(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.ux_max = sortrows(sensor_data.ux_max, 2);

    % remove the reordering data
    sensor_data.ux_max = sensor_data.ux_max(:, 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_max(:, 2) = reorder_index;      
        sensor_data.uy_max = sortrows(sensor_data.uy_max, 2);
        sensor_data.uy_max = sensor_data.uy_max(:, 1);
    end    
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_max(:, 2) = reorder_index;      
        sensor_data.uz_max = sortrows(sensor_data.uz_max, 2);
        sensor_data.uz_max = sensor_data.uz_max(:, 1);
    end      
end

if record.u_min
    % append the reordering data
    sensor_data.ux_min(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.ux_min = sortrows(sensor_data.ux_min, 2);

    % remove the reordering data
    sensor_data.ux_min = sensor_data.ux_min(:, 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_min(:, 2) = reorder_index;      
        sensor_data.uy_min = sortrows(sensor_data.uy_min, 2);
        sensor_data.uy_min = sensor_data.uy_min(:, 1);
    end    
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_min(:, 2) = reorder_index;      
        sensor_data.uz_min = sortrows(sensor_data.uz_min, 2);
        sensor_data.uz_min = sensor_data.uz_min(:, 1);
    end      
end

if record.u_rms
    % append the reordering data
    sensor_data.ux_rms(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.ux_rms = sortrows(sensor_data.ux_rms, 2);

    % remove the reordering data
    sensor_data.ux_rms = sensor_data.ux_rms(:, 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.uy_rms(:, 2) = reorder_index;      
        sensor_data.uy_rms = sortrows(sensor_data.uy_rms, 2);
        sensor_data.uy_rms = sensor_data.uy_rms(:, 1);
    end    
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.uz_rms(:, 2) = reorder_index;      
        sensor_data.uz_rms = sortrows(sensor_data.uz_rms, 2);
        sensor_data.uz_rms = sensor_data.uz_rms(:, 1);
    end      
end

if record.I
    % calculate the position for reordering data
    new_col_pos = length(sensor_data.Ix(1,:)) + 1;        
    
    % append the reordering data
    sensor_data.Ix(:, new_col_pos) = reorder_index;

    % reorder based on the order_index
    sensor_data.Ix = sortrows(sensor_data.Ix, new_col_pos);
    
    % remove the reordering data
    sensor_data.Ix = sensor_data.Ix(:, 1:new_col_pos - 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.Iy(:, new_col_pos) = reorder_index;       
        sensor_data.Iy = sortrows(sensor_data.Iy, new_col_pos);
        sensor_data.Iy = sensor_data.Iy(:, 1:new_col_pos - 1);
    end
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.Iz(:, new_col_pos) = reorder_index; 
        sensor_data.Iz = sortrows(sensor_data.Iz, new_col_pos);
        sensor_data.Iz = sensor_data.Iz(:, 1:new_col_pos - 1);           
    end       
end

if record.I_avg
    % append the reordering data
    sensor_data.Ix_avg(:, 2) = reorder_index;

    % reorder based on the order_index
    sensor_data.Ix_avg = sortrows(sensor_data.Ix_avg, 2);

    % remove the reordering data
    sensor_data.Ix_avg = sensor_data.Ix_avg(:, 1);
    
    % repeat for y-dimension if 2D or 3D
    if kgrid.dim > 1
        sensor_data.Iy_avg(:, 2) = reorder_index;      
        sensor_data.Iy_avg = sortrows(sensor_data.Iy_avg, 2);
        sensor_data.Iy_avg = sensor_data.Iy_avg(:, 1);
    end    
    
    % repeat for z-dimension if 3D
    if kgrid.dim > 2
        sensor_data.Iz_avg(:, 2) = reorder_index;      
        sensor_data.Iz_avg = sortrows(sensor_data.Iz_avg, 2);
        sensor_data.Iz_avg = sensor_data.Iz_avg(:, 1);
    end      
end