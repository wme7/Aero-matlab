function sensor_data = kspaceFirstOrder3DC(varargin)
%KSPACEFIRSTORDER3DC   3D time-domain simulation of wave propagation using C++ code.
%
% DESCRIPTION:
%       kspaceFirstOrder3DC provides an blind interface to the C++ version
%       of kspaceFirstOrder3D (called kspaceFirstOrder3D-OMP). The function
%       works by appending the optional input 'SaveToDisk' to the user
%       inputs and then calling kspaceFirstOrder3D to save the input files
%       to disk. The contents of sensor.record (if set) are parsed as input
%       flags, and the C++ code is run using the system command. The output
%       files are then automatically loaded from disk and returned in the
%       same fashion as kspaceFirstOrder3D. The input and output files are
%       saved to the temporary directory native to the operating system,
%       and are deleted after the function runs.
%
%       This function requires the C++ binary/executable of
%       kspaceFirstOrder3D-OMP to be downloaded from
%       http://www.k-wave.org/download.php and placed in the "binaries"
%       directory of the k-Wave toolbox. Alternatively, the name and
%       location of the binary can be specified using the optional input
%       parameters 'BinaryName' and 'BinariesPath'.
% 
%       Note, not all input options are currently supported, and all
%       display options are ignored (only command line outputs are given).
%       See the k-Wave user manual for more information.
%
%       This function is not recommended for large simulations, as the
%       input variables will reside twice in main memory (once in MATLAB,
%       and once in C++). For large simulations, the C++ code should be
%       called outside of MATLAB. See the k-Wave manual for more
%       information.
%
% USAGE:
%       see kspaceFirstOrder3D
%
% OPTIONAL INPUTS:
%       Optional 'string', value pairs that may be used to modify the
%       default computational settings.
%
%       'BinaryName'    - the name of the binary file (default =
%                         kspaceFirstOrder3D-OMP on linux,
%                         kspaceFirstOrder3D-OMP.exe on windows) 
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
%       author          - Bradley Treeby and Jiri Jaros
%       date            - 3rd February 2012
%       last update     - 28th August 2014
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also kspaceFirstOrder3D, kspaceFirstOrder3DG

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

%#ok<*AGROW>

% extract the optional input arguments 
if nargin > 4
    input_args = varargin(5:end);
else
    input_args = {};
end

% check for a user defined location for the binary
if any(strcmp('BinaryPath', input_args))
   
    % find the location of the optional input in the list
    index = find(strcmp('BinaryPath', input_args));
    
    % read the value of the parameter
    binary_path = input_args{index + 1};
    
    % check for a trailing slash
    if ~strcmp(binary_path(end), filesep)
        binary_path = [binary_path filesep];
    end  
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default path
    binary_path = getkWavePath('binaries');
    
end
    
% check for a user defined name for the binary
if any(strcmp('BinaryName', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('BinaryName', input_args));
    
    % read the value of the parameter
    binary_name = input_args{index + 1}; 
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default name for the binary
    if isunix
        binary_name = 'kspaceFirstOrder3D-OMP';
    else
        binary_name = 'kspaceFirstOrder3D-OMP.exe';
    end
    
end

% check the binaries exist and are in the correct place before doing
% anything else
if ~exist([binary_path, binary_name], 'file')
    
    error(['The binary file ' binary_name ' could not be found in ' binary_path ...
        '. To use the C++ code, the C++ binaries for your operating system must be downloaded from www.k-wave.org/download.php and placed in the binaries folder.']);
    
end

% check for a user defined location for the input and output files 
if any(strcmp('DataPath', input_args))
   
    % find the location of the optional input in the list
    index = find(strcmp('DataPath', input_args));
    
    % read the value of the parameter
    data_path = input_args{index + 1};
    
    % check for a trailing slash
    if ~strcmp(data_path(end), filesep)
        data_path = [data_path filesep];
    end
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set default path
    data_path = tempdir;
    
end
    
% check for a user defined name for the input and output files
if any(strcmp('DataName', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('DataName', input_args));
    
    % read the value of the parameter
    name_prefix = input_args{index + 1};
    input_filename = [name_prefix '_input.h5'];
    output_filename = [name_prefix '_output.h5'];
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set the filename inputs to store data in the default temp directory
    date_string = getDateString;
    input_filename = ['kwave_input_data' date_string '.h5'];
    output_filename = ['kwave_output_data' date_string '.h5'];  
    
end

% add pathname to input and output filenames
input_filename = [data_path, input_filename];
output_filename = [data_path, output_filename];

% check for delete data input
if any(strcmp('DeleteData', input_args))
    
    % find the location of the optional input in the list
    index = find(strcmp('DeleteData', input_args));
    
    % read the value of the parameter
    delete_data = input_args{index + 1};
    
    % check option is true or false
    if ~islogical(delete_data)
        error('Optional input ''DeleteData'' must be Boolean');
    end
    
    % delete the optional input from the input arguments list
    input_args(index:index+1) = [];
    
else
    
    % set data to be deleted
    delete_data = true;
    
end

% set empty options string
options_string = '';

% assign pseudonyms for input structures
kgrid  = varargin{1};
source = varargin{3};
sensor = varargin{4};

% check if the sensor mask is defined as cuboid corners
if isfield(sensor, 'mask') && ( size(sensor.mask, 1) == (2*kgrid.dim) )
    cuboid_corners = true;
else
    cuboid_corners = false;
end

% check if performing time reversal, and replace inputs to explicitly use a
% source with a dirichlet boundary condition
if isfield(sensor, 'time_reversal_boundary_data')
        
    % define a new source structure
    clear source;
    source.p_mask = sensor.mask;
    source.p = flipdim(sensor.time_reversal_boundary_data, 2);
    source.p_mode = 'dirichlet';
    
    % define a new sensor structure
    clear sensor;
    sensor.mask = ones(varargin{1}.Nx, varargin{1}.Ny, varargin{1}.Nz);
    sensor.record = {'p_final'};
    
    % set time reversal flag
    time_rev = true;

else
    
    % set time reversal flag
    time_rev = false;
    
end

% check if sensor.record is given
if isfield(sensor, 'record')
    
    % set the options string to record the required output fields
    if ismember('p', sensor.record)
        options_string = [options_string ' --p_raw'];
    end
    if ismember('p_max', sensor.record)
        options_string = [options_string ' --p_max'];
    end
    if ismember('p_min', sensor.record)
        options_string = [options_string ' --p_min'];
    end    
    if ismember('p_rms', sensor.record)
        options_string = [options_string ' --p_rms'];
    end 
    if ismember('p_max_all', sensor.record)
        options_string = [options_string ' --p_max_all'];
    end
    if ismember('p_min_all', sensor.record)
        options_string = [options_string ' --p_min_all'];
    end     
    if ismember('p_final', sensor.record)
        options_string = [options_string ' --p_final'];
    end    
    if ismember('u', sensor.record)
        options_string = [options_string ' --u_raw'];
    end
    if ismember('u_max', sensor.record)
        options_string = [options_string ' --u_max'];
    end
    if ismember('u_min', sensor.record)
        options_string = [options_string ' --u_min'];
    end    
    if ismember('u_rms', sensor.record)
        options_string = [options_string ' --u_rms'];
    end
    if ismember('u_max_all', sensor.record)
        options_string = [options_string ' --u_max_all'];
    end
    if ismember('u_min_all', sensor.record)
        options_string = [options_string ' --u_min_all'];
    end     
    if ismember('u_final', sensor.record)
        options_string = [options_string ' --u_final'];
    end      
    if ismember('u_non_staggered', sensor.record) || ismember('I_avg', sensor.record) || ismember('I', sensor.record)
        options_string = [options_string ' --u_non_staggered_raw'];
    end
    if (ismember('I_avg', sensor.record) || ismember('I', sensor.record)) && (~ismember('p', sensor.record))
        options_string = [options_string ' --p_raw'];
    end
    
else
    
    % if sensor.record is not given, record the raw time series of p
    options_string = [options_string ' --p_raw'];
    
end
    
% check if sensor.record_start_imdex is given
if isfield(sensor, 'record_start_index')
    options_string = [options_string ' -s ' num2str(sensor.record_start_index)];
end
 
% append the save to disk parameter
input_args = [input_args {'SaveToDisk', input_filename}];

% run the MATLAB code first to generate the input file and save to disk
kspaceFirstOrder3D(varargin{1:2}, source, sensor, input_args{:});
 
% run the simulation in C++ and print outputs to the MATLAB command line
if isunix
    
    % prepend spaces in linux pathnames with \ to allow cd to work
    binary_path = strrep(binary_path, ' ', '\ ');
    
    % clear the library path to prevent a conflict with the FFTW libraries
    % loaded automatically by MATLAB, and run linux binary
    run_string = ['system(''export LD_LIBRARY_PATH=; cd ' binary_path '; ./' binary_name ' -i ' input_filename ' -o ' output_filename options_string ''' ,''-echo'');'];
    eval(run_string);
    
else
    % run Windows binary, placing the input and output filenames in double
    % quotations to avoid problems with spaces
    run_string = ['system(''cd /d ' binary_path ' & ' binary_name ' -i "' input_filename '" -o "' output_filename '" ' options_string ''' ,''-echo'');'];
    eval(run_string);
end

% get the simulation and pml sizes from the output file
Nx = h5read(output_filename, '/Nx');
Ny = h5read(output_filename, '/Ny');
Nz = h5read(output_filename, '/Nz');
pml_x_size = h5read(output_filename, '/pml_x_size');
pml_y_size = h5read(output_filename, '/pml_y_size');
pml_z_size = h5read(output_filename, '/pml_z_size');

% set the default index variables for the _all and _final variables
x1 = 1;
x2 = Nx;
y1 = 1;
y2 = Ny;
z1 = 1;
z2 = Nz; 

% check if the PML is set to be outside the computational grid
pml_inside_index = find(cellfun(@(x) strcmp(x, 'PMLInside'), input_args));
if ~isempty(pml_inside_index)
    
    % if the option is defined, check if PMLInside is false
    if ~input_args{pml_inside_index + 1};

        % if the PML is outside, set the index variables to remove the pml
        % from the _all and _final variables
        x1 = 1 + pml_x_size;
        x2 = Nx - pml_x_size;
        y1 = 1 + pml_y_size;
        y2 = Ny - pml_y_size;
        z1 = 1 + pml_z_size;
        z2 = Nz - pml_z_size;        
        
    end  
end

% load the C++ data back from disk using h5read
if time_rev
    sensor_data = h5read(output_filename, '/p_final');
    sensor_data = sensor_data(x1:x2, y1:y2, z1:z2);
elseif isfield(sensor, 'record')
    
    % load data common to both cuboid corners and binary sensor mask
    if ismember('p_max_all', sensor.record)
        sensor_data.p_max_all = h5read(output_filename, '/p_max_all');
        sensor_data.p_max_all = sensor_data.p_max_all(x1:x2, y1:y2, z1:z2);
    end
    if ismember('p_min_all', sensor.record)
        sensor_data.p_min_all = h5read(output_filename, '/p_min_all');
        sensor_data.p_min_all = sensor_data.p_min_all(x1:x2, y1:y2, z1:z2);
    end     
    if ismember('p_final', sensor.record)
        sensor_data.p_final = h5read(output_filename, '/p_final');
        sensor_data.p_final = sensor_data.p_final(x1:x2, y1:y2, z1:z2);
    end  
    if ismember('u_max_all', sensor.record)
        sensor_data.ux_max_all = h5read(output_filename, '/ux_max_all');
        sensor_data.uy_max_all = h5read(output_filename, '/uy_max_all');
        sensor_data.uz_max_all = h5read(output_filename, '/uz_max_all');
        sensor_data.ux_max_all = sensor_data.ux_max_all(x1:x2, y1:y2, z1:z2);
        sensor_data.uy_max_all = sensor_data.uy_max_all(x1:x2, y1:y2, z1:z2);
        sensor_data.uz_max_all = sensor_data.uz_max_all(x1:x2, y1:y2, z1:z2);
    end
    if ismember('u_min_all', sensor.record)
        sensor_data.ux_min_all = h5read(output_filename, '/ux_min_all');
        sensor_data.uy_min_all = h5read(output_filename, '/uy_min_all');
        sensor_data.uz_min_all = h5read(output_filename, '/uz_min_all');
        sensor_data.ux_min_all = sensor_data.ux_min_all(x1:x2, y1:y2, z1:z2);
        sensor_data.uy_min_all = sensor_data.uy_min_all(x1:x2, y1:y2, z1:z2);
        sensor_data.uz_min_all = sensor_data.uz_min_all(x1:x2, y1:y2, z1:z2);        
    end       
    if ismember('u_final', sensor.record)
        sensor_data.ux_final = h5read(output_filename, '/ux_final');
        sensor_data.uy_final = h5read(output_filename, '/uy_final');
        sensor_data.uz_final = h5read(output_filename, '/uz_final');
        sensor_data.ux_final = sensor_data.ux_final(x1:x2, y1:y2, z1:z2);
        sensor_data.uy_final = sensor_data.uy_final(x1:x2, y1:y2, z1:z2);
        sensor_data.uz_final = sensor_data.uz_final(x1:x2, y1:y2, z1:z2);
    end
    
    % load remaining data depending on how the sensor mask is defined
    if ~cuboid_corners
        
        % load data from a binary sensor mask
        if ismember('p', sensor.record)
            sensor_data.p = h5read(output_filename, '/p');
        end
        if ismember('p_max', sensor.record)
            sensor_data.p_max = h5read(output_filename, '/p_max');
        end
        if ismember('p_min', sensor.record)
            sensor_data.p_min = h5read(output_filename, '/p_min');
        end    
        if ismember('p_rms', sensor.record)
            sensor_data.p_rms = h5read(output_filename, '/p_rms');
        end
        if ismember('u', sensor.record)
            sensor_data.ux = h5read(output_filename, '/ux');
            sensor_data.uy = h5read(output_filename, '/uy');
            sensor_data.uz = h5read(output_filename, '/uz');
        end
        if ismember('u_max', sensor.record)
            sensor_data.ux_max = h5read(output_filename, '/ux_max');
            sensor_data.uy_max = h5read(output_filename, '/uy_max');
            sensor_data.uz_max = h5read(output_filename, '/uz_max');
        end
        if ismember('u_min', sensor.record)
            sensor_data.ux_min = h5read(output_filename, '/ux_min');
            sensor_data.uy_min = h5read(output_filename, '/uy_min');
            sensor_data.uz_min = h5read(output_filename, '/uz_min');
        end    
        if ismember('u_rms', sensor.record)
            sensor_data.ux_rms = h5read(output_filename, '/ux_rms');
            sensor_data.uy_rms = h5read(output_filename, '/uy_rms');
            sensor_data.uz_rms = h5read(output_filename, '/uz_rms');
        end 
        if ismember('u_non_staggered', sensor.record) || ismember('I_avg', sensor.record) || ismember('I', sensor.record)
            sensor_data.ux_non_staggered = h5read(output_filename, '/ux_non_staggered');
            sensor_data.uy_non_staggered = h5read(output_filename, '/uy_non_staggered');
            sensor_data.uz_non_staggered = h5read(output_filename, '/uz_non_staggered');
        end
        
    else
        
        % load data from cuboid corners
        for cuboid_index = 1:size(sensor.mask, 2)
           
            if ismember('p', sensor.record)
                sensor_data(cuboid_index).p = h5read(output_filename, ['/p/' num2str(cuboid_index)]); 
            end
            if ismember('p_max', sensor.record)
                sensor_data(cuboid_index).p_max = h5read(output_filename, ['/p_max/' num2str(cuboid_index)]);
            end
            if ismember('p_min', sensor.record)
                sensor_data(cuboid_index).p_min = h5read(output_filename, ['/p_min/' num2str(cuboid_index)]);
            end    
            if ismember('p_rms', sensor.record)
                sensor_data(cuboid_index).p_rms = h5read(output_filename, ['/p_rms/' num2str(cuboid_index)]);
            end
            if ismember('u', sensor.record)
                sensor_data(cuboid_index).ux = h5read(output_filename, ['/ux/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uy = h5read(output_filename, ['/uy/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uz = h5read(output_filename, ['/uz/' num2str(cuboid_index)]);
            end
            if ismember('u_max', sensor.record)
                sensor_data(cuboid_index).ux_max = h5read(output_filename, ['/ux_max/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uy_max = h5read(output_filename, ['/uy_max/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uz_max = h5read(output_filename, ['/uz_max/' num2str(cuboid_index)]);
            end
            if ismember('u_min', sensor.record)
                sensor_data(cuboid_index).ux_min = h5read(output_filename, ['/ux_min/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uy_min = h5read(output_filename, ['/uy_min/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uz_min = h5read(output_filename, ['/uz_min/' num2str(cuboid_index)]);
            end    
            if ismember('u_rms', sensor.record)
                sensor_data(cuboid_index).ux_rms = h5read(output_filename, ['/ux_rms/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uy_rms = h5read(output_filename, ['/uy_rms/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uz_rms = h5read(output_filename, ['/uz_rms/' num2str(cuboid_index)]);
            end 
            if ismember('u_non_staggered', sensor.record) || ismember('I_avg', sensor.record) || ismember('I', sensor.record)
                sensor_data(cuboid_index).ux_non_staggered = h5read(output_filename, ['/ux_non_staggered/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uy_non_staggered = h5read(output_filename, ['/uy_non_staggered/' num2str(cuboid_index)]);
                sensor_data(cuboid_index).uz_non_staggered = h5read(output_filename, ['/uz_non_staggered/' num2str(cuboid_index)]);
            end            
            
        end
        
    end 
else
    if ~cuboid_corners
        sensor_data.p = h5read(output_filename, '/p');
    else
        for cuboid_index = 1:size(sensor.mask, 2)
            sensor_data(cuboid_index).p = h5read(output_filename, ['/p/' num2str(cuboid_index)]); 
        end
    end
end

% compute the intensity outputs
if isfield(sensor, 'record') && (ismember('I_avg', sensor.record) || ismember('I', sensor.record))
  
    % assign sensor variables needed in subscript
    record.cuboid_corners   = cuboid_corners;
    record.I_avg            = ismember('I_avg', sensor.record);
    record.I                = ismember('I', sensor.record);
    record.p                = ismember('p', sensor.record);
    record.u_non_staggered  = ismember('u_non_staggered', sensor.record);
    
    % assign list of cuboid corners
    if record.cuboid_corners
        record.cuboid_corners_list = sensor.mask;
    end
    
    % set flag that the outputs are from the C++ code (cuboid corners is
    % indexed differently)
    save_intensity_matlab_code = false; %#ok<NASGU>
    
    % run subscript to extract intensity values
    kspaceFirstOrder_saveIntensity;
    
end

% filter the recorded time domain pressure signals using a Gaussian filter
% if defined
if ~time_rev && isfield(sensor, 'frequency_response')
    sensor_data.p = gaussianFilter(sensor_data.p, 1/kgrid.dt, sensor.frequency_response(1), sensor.frequency_response(2));
end

% if sensor.record is not given, assign sensor_data.p to sensor_data
if ~isfield(sensor, 'record') && ~cuboid_corners
    sensor_data = sensor_data.p;
end

% delete the input and output files
if delete_data 
    delete(input_filename);
    delete(output_filename);
end