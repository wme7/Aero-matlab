function mat_sm = smooth(kgrid, mat, restore_max, window_type)
%SMOOTH     Smooth a matrix.
%
% DESCRIPTION:
%       smooth filters an input matrix using an n-dimensional
%       frequency domain window created using getWin. If no window type is
%       specified, a Blackman window is used.
%
% USAGE:
%       mat_sm = smooth(kgrid, mat)
%       mat_sm = smooth(kgrid, mat, restore_max)
%       mat_sm = smooth(kgrid, mat, [], window_type)
%       mat_sm = smooth(kgrid, mat, restore_max, window_type)
%
% INPUTS:
%       kgrid       - k-Wave grid structure returned by makeGrid
%       mat         - spatial distribution to smooth
%
% OPTIONAL INPUTS:
%       restore_max - Boolean controlling whether the maximum value is
%                     restored after smoothing (default = false).
%       window_type - shape of the smoothing window; any valid inputs to
%                     getWin are supported (default = 'Blackman').
%
% OUTPUTS:
%       mat_sm      - smoothed spatial distribution
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 29th April 2009
%       last update - 19th July 2011
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also fft, ifft, fft2, ifft2, fftn, ifftn, makeGrid, getWin, numDim

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

% define literals
DEFAULT_WINDOW_TYPE = 'Blackman';
USE_ROTATION = true;
SYMMETRIC_WINDOW = true;

% check optional inputs
if nargin < 3 || isempty(restore_max)
    restore_max = false;
end

if nargin < 4
    window_type = DEFAULT_WINDOW_TYPE;
end

% check for release B.0.1 inputs of smooth(mat, kgrid, restore_max)
if ~strcmp(class(kgrid), 'kWaveGrid') && ~isstruct(kgrid) %#ok<STISA>
    disp('WARNING: smooth function inputs changed, please update usage');
    mat_sm = smooth(mat, kgrid, restore_max);
    return
end

% extract the number of dimensions and create the filter
switch numDim(mat)
    case 1
        win = getWin(kgrid.Nx, window_type, 'Rotation', USE_ROTATION , 'Symmetric', SYMMETRIC_WINDOW);
    case 2
        win = getWin([kgrid.Nx, kgrid.Ny], window_type, 'Rotation', USE_ROTATION , 'Symmetric', SYMMETRIC_WINDOW);
    case 3
        win = getWin([kgrid.Nx, kgrid.Ny, kgrid.Nz], window_type, 'Rotation', USE_ROTATION , 'Symmetric', SYMMETRIC_WINDOW); 
end

% apply the filter
mat_sm = real(ifftn(fftn(mat).*ifftshift(win)));

% restore magnitude if required
if restore_max
    mat_sm = (max(abs(mat(:)))/max(abs(mat_sm(:))))*mat_sm;
end