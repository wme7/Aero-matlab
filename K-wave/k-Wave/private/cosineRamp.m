function ramp = cosineRamp(length)
%RAMP   create a smoothly varying cosine ramp from 0 to 1.
%
% USAGE:
%       ramp = cosineRamp(length)
%
% ABOUT:
%       author: Bradley E. Treeby
%       date: 2nd December 2009
%       last update: 3rd December 2009
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

ramp = (-cos( (0:(length-1))*pi/length ) + 1)/2;