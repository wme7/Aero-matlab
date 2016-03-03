function sensor_data = kspaceFirstOrder3DG(varargin)
%KSPACEFIRSTORDER3DG   3D time-domain simulation of wave propagation on a GPU using C++ CUDA code.
%
% DESCRIPTION:
%       kspaceFirstOrder3DG provides an blind interface to the native C++
%       CUDA version of kspaceFirstOrder3D (called kspaceFirstOrder3D-CUDA)
%       in the same way as kspaceFirstOrder3DC.
%
%       The function works by appending the optional input 'SaveToDisk' to
%       the user inputs and then calling kspaceFirstOrder3D to save the
%       input files to disk. The contents of sensor.record (if set) are
%       parsed as input flags, and the C++ CUDA code is run using the
%       system command. The output files are then automatically loaded from
%       disk and returned in the same fashion as kspaceFirstOrder3D. The
%       input and output files are saved to the temporary directory native
%       to the operating system, and are deleted after the function runs.
%
%       This function requires the C++ binary/executable of
%       kspaceFirstOrder3D-CUDA to be downloaded from
%       http://www.k-wave.org/download.php and placed in the "binaries"
%       directory of the k-Wave toolbox. 
% 
%       Note, not all input options are currently supported, and all
%       display options are ignored (only command line outputs are given).
%       See the k-Wave user manual for more information.
%
% USAGE:
%       see kspaceFirstOrder3D
%
% OPTIONAL INPUTS:
%       Optional 'string', value pairs that may be used to modify the
%       default computational settings.
%
%       'BinaryName'    - the name of the binary file (default =
%                         kspaceFirstOrder3D-GPU on linux,
%                         kspaceFirstOrder3D-GPU.exe on windows) 
%       'BinaryPath'    - path of the binary file (default = binaries/)
%       'DataName'      - prefix used to generate a custom name for the
%                         input and output data files (this is appended
%                         with _input.h5 and _output.h5) (default =
%                         kwave_<input/output>_data_<date>.h5)
%       'DataPath'      - location of the folder where the input and output
%                         HDF5 files should be stored (default = tempdir)
%       'DeleteData'    - Boolean controlling whether the input and output
%                         HDF5 files should be deleted after running the
%                         simulation (default = true)
%
% ABOUT:
%       author          - Bradley Treeby
%       date            - 30th October 2013
%       last update     - 14th August 2014
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also kspaceFirstOrder3D, kspaceFirstOrder3DC

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

% This function is essentially a wrapper and directly uses the capabilities
% of kspaceFirstOrder3DC by replacing the binary name with the name of the
% GPU binary. 

% check for a custom binary name
if any(strcmp('BinaryName', varargin))
    
    % if the binary name is given, directly pass to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:});
    
else

    % if the binary name is not given, specify to use the GPU binary
    if isunix
        binary_name = 'kspaceFirstOrder3D-CUDA';
    else
        binary_name = 'kspaceFirstOrder3D-CUDA.exe';
    end
    
    % pass this and the original inputs to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:}, 'BinaryName', binary_name);
    
end
