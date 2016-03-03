function c = speedSoundWater(T)
%SPEEDSOUNDWATER   Calculate the speed of sound in water with temperature.
%
% DESCRIPTION:
%       speedSoundWater calculates the speed of sound in distilled water at
%       a given temperature using the 5th order polynomial given by
%       Marczak. 
%
% USAGE:
%       c = speedSoundWater(T)
%
% INPUTS:
%       T   - water temperature [degC]
%
% OUTPUTS:
%       c   - speed of sound [m/s]
%
% ABOUT:
%       author: Bradley E. Treeby
%       date: 11th August 2008
%       reference: Marczak (1997) "Water as a standard in the measurements
%       of speed of sound in liquids," J. Acoust. Soc. Am., 102, 2776-2779
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Coxx

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

p(1) = 2.787860e-9;
p(2) = -1.398845e-6;
p(3) = 3.287156e-4;
p(4) = -5.779136e-2;
p(5) = 5.038813;
p(6) = 1.402385e3;
c = polyval(p, T);