function dim = numDim(x)
%NUMDIM   Return the number of matrix dimensions.
%
% DESCRIPTION:
%       numDim returns the number of dimensions of a matrix x. Unlike
%       the inbuilt ndims, singleton dimensions are not counted, so numDim
%       returns 1 for 1D vectors, 2 for 2D arrays, etc. 
%
% USAGE:
%       dim = numDim(x)
%
% INPUTS:
%       x           - matrix with unknown number of dimensions
%
% OUTPUTS:
%       dim         - number of dimensions
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 29th April 2009
%       last update - 29th February 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also ndims

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

% get the size collapsing any singleton dimensions
sz = size(squeeze(x));

if length(sz) > 2
    dim = length(sz);
elseif sz(1) == 1 || sz(2) == 1
    dim = 1;
else
    dim = 2;
end

% % alternate check for 1D vector in any dimension
% if max(sz) == prod(sz)
%     dim = 1;
% end