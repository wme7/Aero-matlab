% DESCRIPTION:
%       subscript to save a movie file
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 25th November 2010
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

% update command line status
disp('  saving movie file...');

% check if user provided a compression setting
search_found = false;
if ~isempty(movie_args)
    search_index = 1;
    while ~search_found && search_index <= length(movie_args)
        if strcmp(movie_args{search_index}, 'Compression')
            search_found = true;
        end
        search_index = search_index + 1;
    end
end

if search_found
    % use compression setting provided by user
    movie2avi(movie_frames, movie_name, movie_args{:});
else
    % check if on 64 bit, UNIX, MAC or Windows computer
    if sum(strfind(computer, '64'))
        movie2avi(movie_frames, movie_name, 'Compression', MOVIE_COMP_64B, movie_args{:});
    elseif sum(strfind(computer, 'GLNX'))
        movie2avi(movie_frames, movie_name, 'Compression', MOVIE_COMP_LNX, movie_args{:});
    elseif sum(strfind(computer, 'MAC'))
        movie2avi(movie_frames, movie_name, 'Compression', MOVIE_COMP_MAC, movie_args{:});
    else
        movie2avi(movie_frames, movie_name, 'Compression', MOVIE_COMP_WIN, movie_args{:});
    end
end