% DESCRIPTION:
%       subscript to initialise figure window
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 11th February 2014
%       last update - 11th February 2014
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

% create figure window
img = figure;

% resize if using the 2D elastic code
if elastic_code && kgrid.dim == 2
    scaleFig(1.5, 1);
end

% create waitbar
if ~record.time_rev
    if ~elastic_code
        pbar = waitbar(0, 'Computing Pressure Field', 'Visible', 'off');
    else
        pbar = waitbar(0, 'Computing Wavefield', 'Visible', 'off');
    end
else
    pbar = waitbar(0, 'Computing Time Reversed Field', 'Visible', 'off');
end

% shift the waitbar so it doesn't overlap the figure window
posn_pbar = get(pbar, 'OuterPosition');
posn_img = get(img, 'OuterPosition');
posn_pbar(2) = max(min(posn_pbar(2) - posn_pbar(4), posn_img(2) - posn_pbar(4) - 10), 0);
set(pbar, 'OuterPosition', posn_pbar, 'Visible', 'on');