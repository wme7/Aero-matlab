function enforceFields(structure, field_names)
%ENFORCEFIELDS   Checks structure field names for existance.
%
% DESCRIPTION:
%       enforceFields checks a MATLAB structure for the existance of a set
%       of required field name defined by a cell array. If a field name
%       within the cell array is not found within the structure, an error
%       is thrown with the name of the missing field given in the error
%       message. The function is useful for enforcing a set of required
%       fields within a structure.
%
% USAGE:
%       enforceFields(structure, field_names)
%
% INPUTS:
%       structure   - MATLAB structure
%       field_names - cell array of allowable field names
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 13th October 2009
%       last update - 13th October 2009
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also fieldnames

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

names = fieldnames(structure);
for field_index = 1:length(field_names)
    field_found = false;
    for names_index = 1:length(names)
        if strcmp(field_names{field_index}, names{names_index})
            field_found = true;
        end
    end
    if ~field_found
    	error(['The field ' inputname(1) '.' field_names{field_index} ' must be defined']);
    end
end