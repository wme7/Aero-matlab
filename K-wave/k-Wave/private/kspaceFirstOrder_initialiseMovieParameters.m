% DESCRIPTION:
%       subscript to initialise movie parameters
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 21st February 2011
%       last update - 9th February 2012
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

% force getframe compatability with dual monitors
movegui(img);

% define frame index
frame_index = 1;

% calculate the total number of movie frames. The additional two frames
% come from storing the first and last time steps. The last frame is
% already included if plot_freq divides into index_end with no remainder.
if rem(index_end, plot_freq)
    n_frames = floor(index_end / plot_freq) + 2;
else
    n_frames = floor(index_end / plot_freq) + 1;
end

% preallocate movie structure.
movie_frames(1:n_frames) = struct('cdata', [], 'colormap', []);