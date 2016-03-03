function scaleFig(horizontal, vertical)
%SCALEFIG   Resize current figure window.
%
% DESCRIPTION:
%       scaleFig resizes the current figure window by the specified
%       horizontal and vertical multipliers. For example, to increase the
%       width of the current figure by 50%, use the syntax scaleFig(1.5, 1)
%
% USAGE:
%       scaleFig(horizontal, vertical)
%
% INPUTS:
%       horizontal      - size multiplier for the horizontal direction
%       vertical        - size multiplier for the vertical direction
%       
% ABOUT:
%       author          - Bradley Treeby
%       date            - 7th July 2011
%       last update     - 17th April 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also set, get

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

% get the current figure and size
fig_size = get(gcf, 'Position');

% shrink the horizontal direction
fig_size(3) = fig_size(3)*horizontal;

% shrink the vertical direction
fig_size(4) = fig_size(4)*vertical;

% apply size
set(gcf, 'Position', fig_size);

% make sure the figure is on screen
movegui(gca);