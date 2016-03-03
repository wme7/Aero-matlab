function overlayPlot(varargin)
%IMAGEOVERLAY  Overlay two images
%
% DESCRIPTION:
%       overlayPlot overlays two 2D images. The background is displayed
%       using a grayscale map. The foreground is log compressed and
%       thresholded to a particular dynamic range, and overlayed using an
%       alpha value of 0.5.
%
%       Example:
%           x = rand(128);
%           y = peaks(128);
%           overlayPlot(x, y, 30);
%
% USAGE:
%       overlayPlot(bg, fg)
%       overlayPlot(bg, fg, fg_dnr)
%       overlayPlot(x, y, bg, fg, fg_dnr)
%
% INPUTS:
%       x, y        - vectors describing the position of the pixels in the
%                     image equivalent to image(x, y, c)
%       bg          - background image
%       fg          - foreground image
%       fg_dbr      - dynamic range in dB of foreground image
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 17th October 2012
%       last update - 4th June 2013
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

% set the literals
num_colors = 128;
transparency = 0.5;
fg_dnr = 30;

% extract the inputs
if nargin == 2 || nargin == 3
    bg = varargin{1};
    fg = varargin{2};
    if nargin == 3
        fg_dnr = varargin{3};
    end
elseif nargin == 5
    x_vec = varargin{1};
    y_vec = varargin{2};
    bg = varargin{3};
    fg = varargin{4};
    fg_dnr = varargin{5};
else
    error('incorrect number of inputs');
end

% scale the background image from 0 to num_colors
bg = bg - min(bg(:));
bg = round(num_colors*bg/max(bg(:)));

% convert the background image to true color
bg = ind2rgb(bg, gray(num_colors));

% plot the background image
if nargin == 5
    image(x_vec, y_vec, bg);
else
    image(bg);
end

% if a value for the dynamic range is given, log compress and threshold
if fg_dnr
    fg(fg <= 0) = 0;
    fg = 20*log10(fg./max(fg(:)));
    fg = fg + fg_dnr;
    fg(fg <= 0) = 0;
end

% scale the background image from 0 to num_colors
fg = fg - min(fg(:));
fg = round(num_colors*fg/max(fg(:)));

% compute the alpha channel
alpha = transparency*ones(size(fg));
alpha(fg == 0) = 0;

% convert the background image to true color
fg = ind2rgb(fg, jet(num_colors));

% plot the foreground image and set the alpha channel
hold on;
if nargin == 5
    fg_im = image(x_vec, y_vec, fg);
else
    fg_im = image(fg);
end
set(fg_im, 'AlphaData', alpha);