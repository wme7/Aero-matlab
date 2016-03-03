% DESCRIPTION:
%       subscript to scale source terms to the correct units
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 15th February 2012
%       last update - 20th January 2014
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

% get the dimension size
N = kgrid.dim;

% check for non-uniform grid and give error for source terms that haven't
% yet been implemented
if (nonuniform_grid && ( uy_source || uz_source || transducer_source))
    disp('WARNING: source scaling not implemented for non-uniform grids with given source condition');
    return
end

% Scale the input pressure by 1/c^2 (to convert to units of density), then
% by 1/N (to split the input across the split density field). If the
% pressure is injected as a mass source, also scale the pressure by 
% 2*dt*c/dx to account for the time step and convert to units of 
% [kg/(m^3 s)] 
if p_source 
    if strcmp(source.p_mode, 'dirichlet')
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous
            % sound speed 
            source.p = source.p ./ (N*c^2);
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for p_index = 1:length(source.p(:,1))        
                source.p(p_index, :) = source.p(p_index, :) ./ (N*c(p_source_index(p_index))^2);
            end
        end
    else
        if nonuniform_grid
            
            % create empty matrix
            grid_point_sep = zeros(size(kgrid.x));
            
            % compute averaged grid point seperation map, the interior
            % points are calculated using the average distance to all
            % connected grid points (the edge values are not calculated
            % assuming there are no source points in the PML)
            switch kgrid.dim
                case 1
                    grid_point_sep(2:end-1) = kgrid.x_size*(kgrid.xn(3:end, 1) - kgrid.xn(1:end-2, 1))/2;
                case 2
                    grid_point_sep(2:end-1, 2:end-1) = ...
                        kgrid.x_size*(kgrid.xn(3:end, 2:end-1) - kgrid.xn(1:end-2, 2:end-1))/4 + ...
                        kgrid.y_size*(kgrid.yn(2:end-1, 3:end) - kgrid.yn(2:end-1, 1:end-2))/4;
                case 3
                    grid_point_sep(2:end-1, 2:end-1, 2:end-1) = ...
                        kgrid.x_size*(kgrid.xn(3:end, 2:end-1, 2:end-1) - kgrid.xn(1:end-2, 2:end-1, 2:end-1))/6 + ...
                        kgrid.y_size*(kgrid.yn(2:end-1, 3:end, 2:end-1) - kgrid.yn(2:end-1, 1:end-2, 2:end-1))/6 + ...
                        kgrid.z_size*(kgrid.zn(2:end-1, 2:end-1, 3:end) - kgrid.zn(2:end-1, 2:end-1, 1:end-2))/6;
            end
            
            % compute and apply scale parameter
            for p_index = 1:length(source.p(:,1))
                if numel(c) == 1
                    % compute the scale parameter based on the homogeneous sound speed
                    source.p(p_index, :) = source.p(p_index, :) .* (2.*dt./(N*c*grid_point_sep(p_source_index(p_index))));
                else
                    % compute the scale parameter based on the sound speed at that position
                    source.p(p_index, :) = source.p(p_index, :) .* (2.*dt./(N*c(p_source_index(p_index)).*grid_point_sep(p_source_index(p_index))));
                end
            end
            
            % clear unused variables
            clear grid_point_sep;
            
        else
            if numel(c) == 1
                % compute the scale parameter based on the homogeneous
                % sound speed 
                source.p = source.p .* (2*dt./(N*c*kgrid.dx));
            else
                % compute the scale parameter seperately for each source
                % position based on the sound speed at that position
                for p_index = 1:length(source.p(:,1))        
                    source.p(p_index, :) = source.p(p_index, :) .* (2.*dt./(N*c(p_source_index(p_index)).*kgrid.dx));
                end
            end
        end
    end
end

% scale the stress source by 1/N to divide amoungst the split field
% components, and if source.s_mode is not set to 'dirichlet', also scale by 
% 2*dt*c/dx to account for the time step and convert to units of 
% [kg/(m^3 s)] 
if sxx_source
    if strcmp(source.s_mode, 'dirichlet') || p0_source
        source.sxx = source.sxx ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxx = source.sxx .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxx(:,1))        
                source.sxx(s_index, :) = source.sxx(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end
if syy_source
    if strcmp(source.s_mode, 'dirichlet') || p0_source
        source.syy = source.syy ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.syy = source.syy .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.syy(:,1))        
                source.syy(s_index, :) = source.syy(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end
if szz_source
    if strcmp(source.s_mode, 'dirichlet') || p0_source
        source.szz = source.szz ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.szz = source.szz .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.szz(:,1))        
                source.szz(s_index, :) = source.szz(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end
if sxy_source
    if strcmp(source.s_mode, 'dirichlet')
        source.sxy = source.sxy ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxy = source.sxy .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxy(:,1))        
                source.sxy(s_index, :) = source.sxy(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end
if sxz_source
    if strcmp(source.s_mode, 'dirichlet')
        source.sxz = source.sxz ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.sxz = source.sxz .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.sxz(:,1))        
                source.sxz(s_index, :) = source.sxz(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end
if syz_source
    if strcmp(source.s_mode, 'dirichlet')
        source.syz = source.syz ./ N;
    else
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound
            % speed  
            source.syz = source.syz .* (2*dt*c./(N*kgrid.dx));
        else
            % compute the scale parameter seperately for each source
            % position based on the sound speed at that position
            for s_index = 1:length(source.syz(:,1))        
                source.syz(s_index, :) = source.syz(s_index, :) .* (2.*dt.*c(s_source_index(s_index))./(N*kgrid.dx));
            end
        end         
    end
end

% if source.u_mode is not set to 'dirichlet', scale the x-direction
% velocity source terms by 2*dt*c/dx to account for the time step and
% convert to units of [m/s^2] 
if ux_source && ~strcmp(source.u_mode, 'dirichlet')
    
    if nonuniform_grid
    
        % create empty matrix
        grid_point_sep = zeros(size(kgrid.x));
        
        % compute averaged grid point seperation map, the interior
        % points are calculated using the average distance to all
        % connected grid points (the edge values are not calculated
        % assuming there are no source points in the PML)
        grid_point_sep(2:end-1, :, :) = kgrid.x_size*(kgrid.xn(3:end, :, :) - kgrid.xn(1:end-2, :, :))/2;
        
        % compute and apply scale parameter
        for u_index = 1:length(source.ux(:,1))
            if numel(c) == 1
                % compute the scale parameter based on the homogeneous sound speed
                source.ux(u_index, :) = source.ux(u_index, :) .* (2*c*dt./(grid_point_sep(u_source_index(u_index))));
            else
                % compute the scale parameter based on the sound speed at that position
                source.ux(u_index, :) = source.ux(u_index, :) .* (2*c(u_source_index(u_index))*dt./(grid_point_sep(u_source_index(u_index))));
            end
        end
        
        % clear unused variables
        clear grid_point_sep;
    
    else
        
        if numel(c) == 1
            % compute the scale parameter based on the homogeneous sound speed
            source.ux = source.ux .* (2*c*dt./kgrid.dx);
        else
            % compute the scale parameter seperately for each source position
            % based on the sound speed at that position
            for u_index = 1:length(source.ux(:,1))
                source.ux(u_index, :) = source.ux(u_index, :) .* (2.*c(u_source_index(u_index)).*dt./kgrid.dx);
            end
        end
    end
end

% if source.u_mode is not set to 'dirichlet', scale the y-direction
% velocity source terms by 2*dt*c/dy to account for the time step and
% convert to units of [m/s^2] 
if uy_source && ~strcmp(source.u_mode, 'dirichlet')
    if numel(c) == 1
        % compute the scale parameter based on the homogeneous sound speed
        source.uy = source.uy .* (2*c*dt./kgrid.dy);
    else
        % compute the scale parameter seperately for each source position
        % based on the sound speed at that position
        for u_index = 1:length(source.uy(:,1))
            source.uy(u_index, :) = source.uy(u_index, :) .* (2.*c(u_source_index(u_index)).*dt./kgrid.dy);
        end
    end 
end 

% if source.u_mode is not set to 'dirichlet', scale the z-direction
% velocity source terms by 2*dt*c/dz to account for the time step and
% convert to units of [m/s^2]  
if uz_source && ~strcmp(source.u_mode, 'dirichlet') 
    if numel(c) == 1
        % compute the scale parameter based on the homogeneous sound speed
        source.uz = source.uz .* (2*c*dt./kgrid.dz);
    else
        % compute the scale parameter seperately for each source position
        % based on the sound speed at that position
        for u_index = 1:length(source.uz(:,1))        
            source.uz(u_index, :) = source.uz(u_index, :) .* (2.*c(u_source_index(u_index)).*dt./kgrid.dz);
        end
    end
end

% scale the transducer source term by 2*dt*c/dx to account for the time
% step and convert to units of [m/s^2] 
if transducer_source   
    if numel(c) == 1
        transducer_input_signal = transducer_input_signal .* (2*c*dt./kgrid.dx);
    else
        % compute the scale parameter based on the average sound speed at the
        % transducer positions (only one input signal is used to drive the
        % transducer)
        transducer_input_signal = transducer_input_signal.* (2*(mean(c(u_source_index)))*dt./kgrid.dx);
    end
end

% clear subscript variables
clear N;