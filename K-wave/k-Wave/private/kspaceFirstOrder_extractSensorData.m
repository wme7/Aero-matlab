function sensor_data = kspaceFirstOrder_extractSensorData(dim, sensor_data, file_index, sensor_mask_index, record, p, ux_sgx, uy_sgy, uz_sgz)
% DESCRIPTION:
%       Function to extract the required sensor data from the acoustic
%       variables at each time step. This is defined as a function rather
%       than a script to avoid the computational overhead of using scripts
%       to access variables local to another function. For k-Wave < V1.1,
%       this code was included directly in the simulation functions.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 9th July 2013
%       last update - 6th September 2013
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
% GRID STAGGERING
% =========================================================================

% shift the components of the velocity field onto the non-staggered
% grid if required for output
if (...
        record.u_non_staggered || ...
        record.I || ...
        record.I_avg)

    switch dim
        case 1
            ux_shifted = real(ifft(bsxfun(@times, record.x_shift_neg, fft(ux_sgx, [], 1)), [], 1));
        case 2
            ux_shifted = real(ifft(bsxfun(@times, record.x_shift_neg, fft(ux_sgx, [], 1)), [], 1));
            uy_shifted = real(ifft(bsxfun(@times, record.y_shift_neg, fft(uy_sgy, [], 2)), [], 2));
        case 3
            ux_shifted = real(ifft(bsxfun(@times, record.x_shift_neg, fft(ux_sgx, [], 1)), [], 1));
            uy_shifted = real(ifft(bsxfun(@times, record.y_shift_neg, fft(uy_sgy, [], 2)), [], 2));
            uz_shifted = real(ifft(bsxfun(@times, record.z_shift_neg, fft(uz_sgz, [], 3)), [], 3));
    end
end

% =========================================================================
% BINARY SENSOR MASK
% =========================================================================

if record.binary_sensor_mask
            
    % store the time history of the acoustic pressure
    if record.p || record.I || record.I_avg
        if ~record.compute_directivity
            sensor_data.p(:, file_index) = p(sensor_mask_index);
        end
    end

    % store the maximum acoustic pressure
    if record.p_max
        if file_index == 1
            sensor_data.p_max = p(sensor_mask_index);
        else
            sensor_data.p_max = max(sensor_data.p_max, p(sensor_mask_index));
        end
    end

    % store the minimum acoustic pressure
    if record.p_min
        if file_index == 1
            sensor_data.p_min = p(sensor_mask_index);
        else
            sensor_data.p_min = min(sensor_data.p_min, p(sensor_mask_index));
        end                
    end

    % store the rms acoustic pressure
    if record.p_rms
        sensor_data.p_rms = sqrt((sensor_data.p_rms.^2*(file_index - 1) + p(sensor_mask_index).^2)./file_index);
    end
    
    % store the time history of the particle velocity on the staggered grid
    if record.u
        switch dim
            case 1
                sensor_data.ux(:, file_index) = ux_sgx(sensor_mask_index);
            case 2
                sensor_data.ux(:, file_index) = ux_sgx(sensor_mask_index);
                sensor_data.uy(:, file_index) = uy_sgy(sensor_mask_index);                           
            case 3
                sensor_data.ux(:, file_index) = ux_sgx(sensor_mask_index);
                sensor_data.uy(:, file_index) = uy_sgy(sensor_mask_index);            
                sensor_data.uz(:, file_index) = uz_sgz(sensor_mask_index);  
        end
    end
    
    % store the time history of the particle velocity
    if record.u_non_staggered || record.I || record.I_avg
        switch dim
            case 1
                sensor_data.ux_non_staggered(:, file_index) = ux_shifted(sensor_mask_index);
            case 2
                sensor_data.ux_non_staggered(:, file_index) = ux_shifted(sensor_mask_index);
                sensor_data.uy_non_staggered(:, file_index) = uy_shifted(sensor_mask_index);                           
            case 3
                sensor_data.ux_non_staggered(:, file_index) = ux_shifted(sensor_mask_index);
                sensor_data.uy_non_staggered(:, file_index) = uy_shifted(sensor_mask_index);            
                sensor_data.uz_non_staggered(:, file_index) = uz_shifted(sensor_mask_index);  
        end
    end

    % store the maximum particle velocity
    if record.u_max
        if file_index == 1
            switch dim
                case 1
                    sensor_data.ux_max = ux_sgx(sensor_mask_index);
                case 2
                    sensor_data.ux_max = ux_sgx(sensor_mask_index);
                    sensor_data.uy_max = uy_sgy(sensor_mask_index);
                case 3
                    sensor_data.ux_max = ux_sgx(sensor_mask_index);
                    sensor_data.uy_max = uy_sgy(sensor_mask_index);
                    sensor_data.uz_max = uz_sgz(sensor_mask_index);
            end
        else
            switch dim
                case 1
                    sensor_data.ux_max = max(sensor_data.ux_max, ux_sgx(sensor_mask_index));
                case 2
                    sensor_data.ux_max = max(sensor_data.ux_max, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_max = max(sensor_data.uy_max, uy_sgy(sensor_mask_index)); 
                case 3
                    sensor_data.ux_max = max(sensor_data.ux_max, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_max = max(sensor_data.uy_max, uy_sgy(sensor_mask_index));            
                    sensor_data.uz_max = max(sensor_data.uz_max, uz_sgz(sensor_mask_index));   
            end
        end
    end        

    % store the minimum particle velocity
    if record.u_min
        if file_index == 1
            switch dim
                case 1
                    sensor_data.ux_min = ux_sgx(sensor_mask_index);                    
                case 2
                    sensor_data.ux_min = ux_sgx(sensor_mask_index);
                    sensor_data.uy_min = uy_sgy(sensor_mask_index);                    
                case 3
                    sensor_data.ux_min = ux_sgx(sensor_mask_index);
                    sensor_data.uy_min = uy_sgy(sensor_mask_index);
                    sensor_data.uz_min = uz_sgz(sensor_mask_index);
            end
        else
            switch dim
                case 1
                    sensor_data.ux_min = min(sensor_data.ux_min, ux_sgx(sensor_mask_index));                      
                case 2
                    sensor_data.ux_min = min(sensor_data.ux_min, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_min = min(sensor_data.uy_min, uy_sgy(sensor_mask_index));                      
                case 3
                    sensor_data.ux_min = min(sensor_data.ux_min, ux_sgx(sensor_mask_index));            
                    sensor_data.uy_min = min(sensor_data.uy_min, uy_sgy(sensor_mask_index));            
                    sensor_data.uz_min = min(sensor_data.uz_min, uz_sgz(sensor_mask_index)); 
            end                       
        end
    end             

    % store the rms particle velocity
    if record.u_rms
        switch dim
            case 1
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + ux_sgx(sensor_mask_index).^2)./file_index);
            case 2
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + ux_sgx(sensor_mask_index).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + uy_sgy(sensor_mask_index).^2)./file_index);
            case 3
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + ux_sgx(sensor_mask_index).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + uy_sgy(sensor_mask_index).^2)./file_index);
                sensor_data.uz_rms(:) = sqrt((sensor_data.uz_rms(:).^2*(file_index - 1) + uz_sgz(sensor_mask_index).^2)./file_index);
        end        
    end   

% =========================================================================
% CARTESIAN SENSOR MASK
% =========================================================================
    
% extract data from specified Cartesian coordinates using interpolation
% (record.record.tri and record.record.bc are the Delaunay triangulation
% and Barycentric coordinates returned by gridDataFast3D)
else

    % store the time history of the acoustic pressure
    if record.p || record.I || record.I_avg
        if dim == 1
            sensor_data.p(:, file_index) = interp1(record.grid_x, p, record.sensor_x);
        else
            sensor_data.p(:, file_index) = sum(p(record.tri) .* record.bc, 2);
        end
    end

    % store the maximum acoustic pressure
    if record.p_max
        if dim == 1
            if file_index == 1
                sensor_data.p_max = interp1(record.grid_x, p, record.sensor_x);
            else
                sensor_data.p_max = max(sensor_data.p_max, interp1(record.grid_x, p, record.sensor_x));
            end            
        else
            if file_index == 1
                sensor_data.p_max = sum(p(record.tri) .* record.bc, 2);
            else
                sensor_data.p_max = max(sensor_data.p_max, sum(p(record.tri) .* record.bc, 2));
            end
        end
    end        

    % store the minimum acoustic pressure
    if record.p_min
        if dim == 1
            if file_index == 1
                sensor_data.p_min = interp1(record.grid_x, p, record.sensor_x);
            else
                sensor_data.p_min = min(sensor_data.p_min, interp1(record.grid_x, p, record.sensor_x));
            end            
        else
            if file_index == 1
                sensor_data.p_min = sum(p(record.tri) .* record.bc, 2);
            else
                sensor_data.p_min = min(sensor_data.p_min, sum(p(record.tri) .* record.bc, 2));
            end
        end
    end            

    % store the rms acoustic pressure
    if record.p_rms
        if dim == 1
            sensor_data.p_rms = sqrt((sensor_data.p_rms.^2*(file_index - 1) + (interp1(record.grid_x, p, record.sensor_x)).^2)./file_index);
        else
            sensor_data.p_rms(:) = sqrt((sensor_data.p_rms(:).^2*(file_index - 1) + (sum(p(record.tri) .* record.bc, 2)).^2)./file_index);
        end
    end

    % store the time history of the particle velocity on the staggered grid
    if record.u
        switch dim
            case 1
                sensor_data.ux(:, file_index) = interp1(record.grid_x, ux_sgx, record.sensor_x);
            case 2
                sensor_data.ux(:, file_index) = sum(ux_sgx(record.tri) .* record.bc, 2);
                sensor_data.uy(:, file_index) = sum(uy_sgy(record.tri) .* record.bc, 2);                
            case 3
                sensor_data.ux(:, file_index) = sum(ux_sgx(record.tri) .* record.bc, 2);
                sensor_data.uy(:, file_index) = sum(uy_sgy(record.tri) .* record.bc, 2);
                sensor_data.uz(:, file_index) = sum(uz_sgz(record.tri) .* record.bc, 2);
        end
    end

    % store the time history of the particle velocity
    if record.u_non_staggered || record.I || record.I_avg
        switch dim
            case 1
                sensor_data.ux_non_staggered(:, file_index) = interp1(record.grid_x, ux_shifted, record.sensor_x);
            case 2
                sensor_data.ux_non_staggered(:, file_index) = sum(ux_shifted(record.tri) .* record.bc, 2);
                sensor_data.uy_non_staggered(:, file_index) = sum(uy_shifted(record.tri) .* record.bc, 2);                
            case 3
                sensor_data.ux_non_staggered(:, file_index) = sum(ux_shifted(record.tri) .* record.bc, 2);
                sensor_data.uy_non_staggered(:, file_index) = sum(uy_shifted(record.tri) .* record.bc, 2);
                sensor_data.uz_non_staggered(:, file_index) = sum(uz_shifted(record.tri) .* record.bc, 2);
        end
    end
    
    % store the maximum particle velocity
    if record.u_max
        if file_index == 1
            switch dim
                case 1
                    sensor_data.ux_max = interp1(record.grid_x, ux_sgx, record.sensor_x);
                case 2
                    sensor_data.ux_max = sum(ux_sgx(record.tri) .* record.bc, 2);
                    sensor_data.uy_max = sum(uy_sgy(record.tri) .* record.bc, 2);                    
                case 3
                    sensor_data.ux_max = sum(ux_sgx(record.tri) .* record.bc, 2);
                    sensor_data.uy_max = sum(uy_sgy(record.tri) .* record.bc, 2);
                    sensor_data.uz_max = sum(uz_sgz(record.tri) .* record.bc, 2);
            end
        else
            switch dim
                case 1
                    sensor_data.ux_max = max(sensor_data.ux_max, interp1(record.grid_x, ux_sgx, record.sensor_x));
                case 2
                    sensor_data.ux_max = max(sensor_data.ux_max, sum(ux_sgx(record.tri) .* record.bc, 2));
                    sensor_data.uy_max = max(sensor_data.uy_max, sum(uy_sgy(record.tri) .* record.bc, 2));
                case 3
                    sensor_data.ux_max = max(sensor_data.ux_max, sum(ux_sgx(record.tri) .* record.bc, 2));
                    sensor_data.uy_max = max(sensor_data.uy_max, sum(uy_sgy(record.tri) .* record.bc, 2));
                    sensor_data.uz_max = max(sensor_data.uz_max, sum(uz_sgz(record.tri) .* record.bc, 2));
            end
        end
    end   

    % store the minimum particle velocity
    if record.u_min
        if file_index == 1
            switch dim
                case 1
                    sensor_data.ux_min = interp1(record.grid_x, ux_sgx, record.sensor_x);
                case 2
                    sensor_data.ux_min = sum(ux_sgx(record.tri) .* record.bc, 2);
                    sensor_data.uy_min = sum(uy_sgy(record.tri) .* record.bc, 2);                    
                case 3
                    sensor_data.ux_min = sum(ux_sgx(record.tri) .* record.bc, 2);
                    sensor_data.uy_min = sum(uy_sgy(record.tri) .* record.bc, 2);
                    sensor_data.uz_min = sum(uz_sgz(record.tri) .* record.bc, 2);
            end
        else
            switch dim
                case 1
                    sensor_data.ux_min = min(sensor_data.ux_min, interp1(record.grid_x, ux_sgx, record.sensor_x));
                case 2
                    sensor_data.ux_min = min(sensor_data.ux_min, sum(ux_sgx(record.tri) .* record.bc, 2));
                    sensor_data.uy_min = min(sensor_data.uy_min, sum(uy_sgy(record.tri) .* record.bc, 2));
                case 3
                    sensor_data.ux_min = min(sensor_data.ux_min, sum(ux_sgx(record.tri) .* record.bc, 2));
                    sensor_data.uy_min = min(sensor_data.uy_min, sum(uy_sgy(record.tri) .* record.bc, 2));
                    sensor_data.uz_min = min(sensor_data.uz_min, sum(uz_sgz(record.tri) .* record.bc, 2));
            end
        end
    end             

    % store the rms particle velocity
    if record.u_rms
        switch dim
            case 1
                sensor_data.ux_rms = sqrt((sensor_data.ux_rms.^2*(file_index - 1) + (interp1(record.grid_x, ux_sgx, record.sensor_x)).^2)./file_index);
            case 2
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + (sum(ux_sgx(record.tri) .* record.bc, 2)).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + (sum(uy_sgy(record.tri) .* record.bc, 2)).^2)./file_index);
            case 3
                sensor_data.ux_rms(:) = sqrt((sensor_data.ux_rms(:).^2*(file_index - 1) + (sum(ux_sgx(record.tri) .* record.bc, 2)).^2)./file_index);
                sensor_data.uy_rms(:) = sqrt((sensor_data.uy_rms(:).^2*(file_index - 1) + (sum(uy_sgy(record.tri) .* record.bc, 2)).^2)./file_index);
                sensor_data.uz_rms(:) = sqrt((sensor_data.uz_rms(:).^2*(file_index - 1) + (sum(uz_sgz(record.tri) .* record.bc, 2)).^2)./file_index);
        end
    end        
end

% =========================================================================
% RECORDED VARIABLES OVER ENTIRE GRID
% =========================================================================

% store the maximum acoustic pressure over all the grid elements
if record.p_max_all
    switch dim
        case 1
            if file_index == 1
                sensor_data.p_max_all = p(record.x1_inside:record.x2_inside);
            else
                sensor_data.p_max_all = max(sensor_data.p_max_all, p(record.x1_inside:record.x2_inside));
            end
        case 2
            if file_index == 1
                sensor_data.p_max_all = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
            else
                sensor_data.p_max_all = max(sensor_data.p_max_all, p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));
            end
        case 3
            if file_index == 1
                sensor_data.p_max_all = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
            else
                sensor_data.p_max_all = max(sensor_data.p_max_all, p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));
            end
    end
end

% store the minimum acoustic pressure over all the grid elements
if record.p_min_all
    switch dim
        case 1
            if file_index == 1
                sensor_data.p_min_all = p(record.x1_inside:record.x2_inside);
            else
                sensor_data.p_min_all = min(sensor_data.p_min_all, p(record.x1_inside:record.x2_inside));
            end              
        case 2
            if file_index == 1
                sensor_data.p_min_all = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
            else
                sensor_data.p_min_all = min(sensor_data.p_min_all, p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));
            end  
        case 3
            if file_index == 1
                sensor_data.p_min_all = p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
            else
                sensor_data.p_min_all = min(sensor_data.p_min_all, p(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));
            end 
    end
end              

% store the maximum particle velocity over all the grid elements
if record.u_max_all
    switch dim
        case 1
            if file_index == 1
                sensor_data.ux_max_all = ux_sgx(record.x1_inside:record.x2_inside);
            else
                sensor_data.ux_max_all = max(sensor_data.ux_max_all, ux_sgx(record.x1_inside:record.x2_inside));            
            end               
        case 2
            if file_index == 1
                sensor_data.ux_max_all = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
                sensor_data.uy_max_all = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
            else
                sensor_data.ux_max_all = max(sensor_data.ux_max_all, ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));            
                sensor_data.uy_max_all = max(sensor_data.uy_max_all, uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));            
            end            
        case 3
            if file_index == 1
                sensor_data.ux_max_all = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
                sensor_data.uy_max_all = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
                sensor_data.uz_max_all = uz_sgz(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
            else
                sensor_data.ux_max_all = max(sensor_data.ux_max_all, ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
                sensor_data.uy_max_all = max(sensor_data.uy_max_all, uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
                sensor_data.uz_max_all = max(sensor_data.uz_max_all, uz_sgz(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
            end
    end
end        

% store the minimum particle velocity over all the grid elements
if record.u_min_all
    switch dim
        case 1
            if file_index == 1
                sensor_data.ux_min_all = ux_sgx(record.x1_inside:record.x2_inside);
            else
                sensor_data.ux_min_all = min(sensor_data.ux_min_all, ux_sgx(record.x1_inside:record.x2_inside));            
            end            
        case 2 
            if file_index == 1
                sensor_data.ux_min_all = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
                sensor_data.uy_min_all = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside);
            else
                sensor_data.ux_min_all = min(sensor_data.ux_min_all, ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));            
                sensor_data.uy_min_all = min(sensor_data.uy_min_all, uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside));                     
            end
        case 3
            if file_index == 1
                sensor_data.ux_min_all = ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
                sensor_data.uy_min_all = uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
                sensor_data.uz_min_all = uz_sgz(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside);
            else
                sensor_data.ux_min_all = min(sensor_data.ux_min_all, ux_sgx(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
                sensor_data.uy_min_all = min(sensor_data.uy_min_all, uy_sgy(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
                sensor_data.uz_min_all = min(sensor_data.uz_min_all, uz_sgz(record.x1_inside:record.x2_inside, record.y1_inside:record.y2_inside, record.z1_inside:record.z2_inside));            
            end
    end
end