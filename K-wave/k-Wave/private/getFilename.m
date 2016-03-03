function [filename, pathname] = getFilename(varargin)
% DESCRIPTION:
%       Function to automate collecting a single filename input using the
%       built in uigetfile gui.
%
% USAGE:
%       [filename pathname] = getFilename(start_dir)
%       [filename pathname] = getFilename()
%
% INPUTS:
%       start_dir   - the directory to start in
%
% OUTPUTS:
%       filename / pathname of selected file
%
% ABOUT:
%       author          - Bradley Treeby
%       date            - 5th November 2008
%       last update     - 30th June 2009
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

% get the filename
if nargin == 1
    start_dir = varargin{1};
    [filename, pathname, filterindex] = uigetfile({'*.jpg; *.gif; *.bmp; *.png', 'Image Files (*.jpg, *.gif, *.bmp, *.png)'; '*.*',...
        'All Files (*.*)'}, 'Select data file', 'MultiSelect', 'off', start_dir);
else
    [filename, pathname, filterindex] = uigetfile({'*.jpg; *.gif; *.bmp; *.png', 'Image Files (*.jpg, *.gif, *.bmp, *.png)'; '*.*',...
        'All Files (*.*)'}, 'Select data file', 'MultiSelect', 'off');
end

% check and end with an error if no file has been selected
if filterindex == 0
    error('Please select a data file');
end