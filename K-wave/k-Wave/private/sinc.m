function y = sinc(x)
%SINC   Mathematical sinc function.     
%
% USAGE:
%       y = sinc(x)
%
% ABOUT:
%       author: Bradley E. Treeby
%       date: 14th January 2009
%       last update: 15th September 2010
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

% warning('off', 'MATLAB:divideByZero');
% y = sin(x)./x;
% y(x == 0) = 1;

% new implementation to avoid changing warning status
zero_vals = (x == 0);
y = sin(x + pi*zero_vals) ./ (x + zero_vals) + zero_vals;