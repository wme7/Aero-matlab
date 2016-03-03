% DESCRIPTION:
%       subscript to set the literals and defaults used in the fluid and
%       elastic simulation codes
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 11th February 2014
%       last update - 25th August 2014
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

%#ok<*NASGU>

% number of input variables required to run the simulation codes
NUM_REQ_INPUT_VARIABLES             = 4;  

% set defaults for literals that can be changed using optional inputs
% (these are used in all codes)
CARTESIAN_INTERP_DEF                = 'linear';
CREATE_LOG_DEF                      = false;
DATA_CAST_DEF                       = 'off';
DATA_RECAST_DEF                     = false;
DISPLAY_MASK_DEF                    = 'default';
LOG_SCALE_DEF                       = false;
LOG_SCALE_COMPRESSION_FACTOR_DEF    = 0.02;
MOVIE_ARGS_DEF                      = {};
MOVIE_NAME_DEF                      = [getDateString '-' MFILE];
PLOT_FREQ_DEF                       = 10;
PLOT_LAYOUT_DEF                     = false;
PLOT_SIM_DEF                        = true;
PLOT_PML_DEF                        = true;
PML_ALPHA_DEF                       = 2;
PML_INSIDE_DEF                      = true;
RECORD_MOVIE_DEF                    = false;
SCALE_SOURCE_TERMS_DEF              = true;
SMOOTH_P0_DEF                       = true;
SMOOTH_C0_DEF                       = false;
SMOOTH_RHO0_DEF                     = false;
SOURCE_S_MODE_DEF                   = 'additive';
SOURCE_P_MODE_DEF                   = 'additive';
SOURCE_U_MODE_DEF                   = 'additive';
USE_KSPACE_DEF                      = true;
USE_SG_DEF                          = true;

% set defaults for literals that can be changed using optional inputs
% (these are only used in 1D)
USE_FINITE_DIFFERENCE_DEF           = false;

% set defaults for literals that can be changed using optional inputs
% (these are only used in 2D)
MESH_PLOT_DEF                       = false;
MOVIE_TYPE_DEF                      = 'frame';
FORCE_TSEARCH                       = false;
DIRECTIVITY_PATTERN_DEF             = 'pressure';
DIRECTIVITY_SIZE_SCALE_FACTOR_DEF   = 10;

% set defaults for literals that can be changed using optional inputs
% (these are only used in 3D)
SAVE_TO_DISK_DEF                    = false;
SAVE_TO_DISK_FILENAME_DEF           = 'kwave_input_data.h5';
SAVE_TO_DISK_EXIT_DEF               = true;
STREAM_TO_DISK_DEF                  = false;
STREAM_TO_DISK_STEPS_DEF            = 200;
STREAM_TO_DISK_FILENAME             = 'temp_sensor_data.bin';

% set default movie compression
MOVIE_COMP_WIN                      = 'Cinepak';
MOVIE_COMP_MAC                      = 'None';
MOVIE_COMP_LNX                      = 'None';
MOVIE_COMP_64B                      = 'None';

% set additional literals that can't be changed using optional inputs
COLOR_MAP                           = getColorMap;
ESTIMATE_SIM_TIME_STEPS             = 50;
LOG_NAME                            = ['k-Wave-Log-' getDateString];

% set additional literals that vary depending on the dimension of the
% simulation code
switch kgrid.dim
    case 1
        PLOT_SCALE_DEF              = [-1.1 1.1];
        PML_SIZE_DEF                = 20;
        PLOT_SCALE_WARNING          = 5;
    case 2
        PLOT_SCALE_DEF              = [-1, 1];
        PML_SIZE_DEF                = 20;
        PLOT_SCALE_WARNING          = 10;
    case 3
        PLOT_SCALE_DEF              = [-1, 1];
        PML_SIZE_DEF                = 10;
        PLOT_SCALE_WARNING          = 20;
end

% set the default CFL values used if kgrid.t_array is set to 'auto'
KSPACE_CFL                          = 0.3;
PSTD_CFL                            = 0.1;

% set additional literals unique to the elastic code
MULTI_AXIAL_PML_RATIO_DEF = 0.1;