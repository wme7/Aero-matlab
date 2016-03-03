function circle = makeCircle(Nx, Ny, cx, cy, radius, arc_angle, plot_circle)
%MAKECIRCLE     Create a binary map of a circle within a 2D grid.
%
% DESCRIPTION:
%       makeCircle creates a binary map of a circle or arc (using the
%       midpoint circle algorithm) within a two-dimensional grid (the
%       circle position is denoted by 1's in the matrix with 0's
%       elsewhere). A single grid point is taken as the circle centre thus
%       the total diameter will always be an odd number of grid points.
%
% USAGE:
%       circle = makeCircle(Nx, Ny, cx, cy, radius)
%       circle = makeCircle(Nx, Ny, cx, cy, radius, arc_angle)
%       circle = makeCircle(Nx, Ny, cx, cy, radius, arc_angle, plot_circle)
%
% INPUTS:
%       Nx, Ny          - size of the 2D grid [grid points]
%       cx, cy          - centre of the circle [grid points], if set
%                         to 0, the centre of the grid is used
%       radius          - circle radius [grid points]
%
% OPTIONAL INPUTS:
%       arc_angle       - arc angle for incomplete circle [radians]
%                         (default = 2*pi)
%       plot_circle     - Boolean controlling whether the circle is plotted
%                         using imagesc (default = false)
%
% OUTPUTS:
%       circle          - 2D binary map of a circle
%
% ABOUT:
%       author          - Bradley Treeby
%       date            - 1st May 2009
%       last update     - 20th December 2011
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also makeCartCircle, makeDisc

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

% check for plot_circle input
if nargin < 7
    plot_circle = false;
end

% check for arc_angle input
if nargin < 6
    arc_angle = 2*pi;
elseif arc_angle > 2*pi
    arc_angle = 2*pi;
elseif arc_angle < 0
    arc_angle = 0;
end

% force integer values
Nx = round(Nx);
Ny = round(Ny);
cx = round(cx);
cy = round(cy);
radius = round(radius);

% check for zero values
if cx == 0
    cx = floor(Nx/2) + 1;
end
if cy == 0
    cy = floor(Ny/2) + 1;
end

% check the inputs
if cx < 1 || cx > Nx || cy < 1 || cy > Ny
    error('The center of the circle must be within the grid');
end

% define literals
MAGNITUDE = 1;

% create empty matrix
circle = zeros(Nx, Ny);

% initialise loop variables
x = 0;
y = radius;
d = 1 - radius;

% draw the first cardinal point
try 
    circle(cx, cy - y) = MAGNITUDE;
catch
    error('The circle must fit within the grid');
end

% draw the remaining cardinal points
py = [cx, cx+y, cx-y];
px = [cy+y, cy, cy];
for point_index = 1:length(py)
    
    % check whether the point is within the arc made by arc_angle
    if (atan2(py(point_index) - cx, px(point_index) - cy) + pi) <= arc_angle
        circle(py(point_index), px(point_index)) = MAGNITUDE;
    end
end

% loop through the remaining points using the midpoint circle algorithm
while ( x < y - 1 )
    
    x = x + 1;
    if ( d < 0 ) 
        d = d + x + x + 1;
    else 
        y = y - 1;
        a = x - y + 1;
        d = d + a + a;
    end
    
    % setup point indices
    py = [x+cx, y+cx, y+cx, x+cx, -x+cx, -y+cx, -y+cx, -x+cx];
    px = [y+cy, x+cy, -x+cy, -y+cy, -y+cy, -x+cy, x+cy, y+cy];
    
    % loop through each point
    for point_index = 1:length(py)
        
        % check whether the point is within the arc made by arc_angle
        if (atan2(py(point_index) - cx, px(point_index) - cy) + pi) <= arc_angle
            circle(py(point_index), px(point_index)) = MAGNITUDE;
        end
    end
end

% create the figure
if plot_circle
    figure;
    imagesc(circle, [-1 1]);
    colormap(getColorMap);
    axis image;
    xlabel('y-position [grid points]');
    ylabel('x-position [grid points]');
end