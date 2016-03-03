% DESCRIPTION:
%       subscript to cast output variables back to double precision
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 1st September 2012
%       last update - 27th August 2014
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
disp('  recasting variables to double...');

% set recast string to cast variables back to double (the parallel
% computing toolbox requires the gather command to be used)
if strncmp(data_cast, 'gpuArray', 8)
    recast_pre = 'double(gather(';
    recast_post = '));';
else
    recast_pre = 'double(';
    recast_post = ');';
end

% time history of the acoustic pressure
if record.p
    eval(['sensor_data.p = ' recast_pre 'sensor_data.p' recast_post]);
end

% maximum pressure
if record.p_max
    eval(['sensor_data.p_max = ' recast_pre 'sensor_data.p_max' recast_post]);
end 

% minimum pressure
if record.p_min
    eval(['sensor_data.p_min = ' recast_pre 'sensor_data.p_min' recast_post]);
end 

% rms pressure
if record.p_rms
    eval(['sensor_data.p_rms = ' recast_pre 'sensor_data.p_rms' recast_post]);
end

% final acoustic pressure over all grid points
if record.p_final
    eval(['sensor_data.p_final = ' recast_pre 'sensor_data.p_final' recast_post]);
end

% maximum pressure over all grid points
if record.p_max_all
    eval(['sensor_data.p_max_all = ' recast_pre 'sensor_data.p_max_all' recast_post]);
end 

% minimum pressure over all grid points
if record.p_min_all
    eval(['sensor_data.p_min_all = ' recast_pre 'sensor_data.p_min_all' recast_post]);
end 

% time history of the acoustic particle velocity
if record.u
    eval(['sensor_data.ux = ' recast_pre 'sensor_data.ux' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy = ' recast_pre 'sensor_data.uy' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz = ' recast_pre 'sensor_data.uz' recast_post]);
    end
end

% time history of the acoustic particle velocity
if record.u
    eval(['sensor_data.ux = ' recast_pre 'sensor_data.ux' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy = ' recast_pre 'sensor_data.uy' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz = ' recast_pre 'sensor_data.uz' recast_post]);
    end
end

% time history of the acoustic particle velocity on non-staggered grid
if record.u_non_staggered
    eval(['sensor_data.ux_non_staggered = ' recast_pre 'sensor_data.ux_non_staggered' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_non_staggered = ' recast_pre 'sensor_data.uy_non_staggered' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_non_staggered = ' recast_pre 'sensor_data.uz_non_staggered' recast_post]);
    end
end

% maximum particle velocity
if record.u_max
    eval(['sensor_data.ux_max = ' recast_pre 'sensor_data.ux_max' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_max = ' recast_pre 'sensor_data.uy_max' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_max = ' recast_pre 'sensor_data.uz_max' recast_post]);
    end
end

% minimum particle velocity
if record.u_min
    eval(['sensor_data.ux_min = ' recast_pre 'sensor_data.ux_min' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_min = ' recast_pre 'sensor_data.uy_min' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_min = ' recast_pre 'sensor_data.uz_min' recast_post]);
    end
end

% rms particle velocity
if record.u_rms
    eval(['sensor_data.ux_rms = ' recast_pre 'sensor_data.ux_rms' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_rms = ' recast_pre 'sensor_data.uy_rms' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_rms = ' recast_pre 'sensor_data.uz_rms' recast_post]);
    end
end

% final particle velocity everywhere within medium
if record.u_final
    eval(['sensor_data.ux_final = ' recast_pre 'sensor_data.ux_final' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_final = ' recast_pre 'sensor_data.uy_final' recast_post]);
    end
    if kgrid.dim > 2
        eval(['sensor_data.uz_final = ' recast_pre 'sensor_data.uz_final' recast_post]);
    end
end

% maximum particle velocity over all grid points
if record.u_max_all
    eval(['sensor_data.ux_max_all = ' recast_pre 'sensor_data.ux_max_all' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_max_all = ' recast_pre 'sensor_data.uy_max_all' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_max_all = ' recast_pre 'sensor_data.uz_max_all' recast_post]);
    end
end

% minimum particle velocity over all grid points
if record.u_min_all
    eval(['sensor_data.ux_min_all = ' recast_pre 'sensor_data.ux_min_all' recast_post]);
    if kgrid.dim > 1
        eval(['sensor_data.uy_min_all = ' recast_pre 'sensor_data.uy_min_all' recast_post]);
    end
    if kgrid.dim > 2    
        eval(['sensor_data.uz_min_all = ' recast_pre 'sensor_data.uz_min_all' recast_post]);
    end
end

% object of the kWaveTransducer class is being used as a sensor
if record.transducer_sensor
    eval(['sensor_data.transducer = ' recast_pre 'sensor_data.transducer' recast_post]);
end

% % time history of the instantaneus acoustic intensity
% if record.I
%     eval(['sensor_data.Ix = ' recast_pre 'sensor_data.Ix' recast_post]);
%     if kgrid.dim > 1
%         eval(['sensor_data.Iy = ' recast_pre 'sensor_data.Iy' recast_post]);
%     end
%     if kgrid.dim > 2
%         eval(['sensor_data.Iz = ' recast_pre 'sensor_data.Iz' recast_post]);   
%     end
% end
% 
% % average acoustic intensity
% if record.I_avg
%     eval(['sensor_data.Ix_avg = ' recast_pre 'sensor_data.Ix_avg' recast_post]);
%     if kgrid.dim > 1
%         eval(['sensor_data.Iy_avg = ' recast_pre 'sensor_data.Iy_avg' recast_post]);
%     end
%     if kgrid.dim > 2
%         eval(['sensor_data.Iz_avg = ' recast_pre 'sensor_data.Iz_avg' recast_post]);
%     end
% end