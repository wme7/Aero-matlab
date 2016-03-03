% DESCRIPTION:
%       subscript to retract transducer grid size
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 26th September 2012
%       last update - 26th September 2012
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

% resize the transducer object if the grid has been expanded
if ~PML_inside

    % set the size to trim from the grid
    retract_size = [PML_x_size, PML_y_size, PML_z_size]; 
    
    % check if the sensor is a transducer
    if isa(sensor, 'kWaveTransducer')
        
        % retract the transducer mask
        sensor.retract_grid(retract_size);
        
    end
        
    % check if the source is a transducer, and if so, and different
    % transducer to the sensor 
    if isa(source, 'kWaveTransducer') && ~(isa(sensor, 'kWaveTransducer') && isequal(sensor, source))
        
        % retract the transducer mask
        source.retract_grid(retract_size);
        
    end
end