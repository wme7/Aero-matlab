function writeMatrix(filename, matrix, matrix_name)
%WRITEMATRIX   Write MATLAB matrix to a k-Wave HDF5 file.
%
% DESCRIPTION:
%       writeMatrix writes a MATLAB matrix to an existing HDF5 file. If the
%       HDF5 file specified by filename doesn't exist, it is automatically
%       created. Matrix attributes required by the k-Wave C++ code are
%       automatically added. All variables are stored as 3D matrices (with
%       the size of unused dimensions set to 1), and complex matrices are
%       reorganised to C++ format before storing.
%
% USAGE:
%       writeMatrix(filename, matrix, matrix_name)
%
% INPUTS:
%       filename    - name of HDF5 file to write matrix to
%       matrix      - data to write
%       matrix_name - string containing the name of the matrix
%
% ABOUT:
%       author      - Bradley Treeby and Jiri Jaros
%       date        - 9th August 2012
%       last update - 26th November 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also h5write, h5compare, writeAttributes, writeFlags, writeGrid

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

% get literals
getH5Literals;

% get the size of the input matrix
[Nx, Ny, Nz] = size(matrix);
dims = numDim(matrix);

% check size of matrix and set chunk size and compression level
switch dims
    
    case 3
    
        % set chunk size to Nx*Ny
        chunk_size = [Nx, Ny, 1];

        % set compression level
        compression_level = COMPRESSION_LEVEL;
    
    case 2
        
        % set chunk size to Nx
        chunk_size = [Nx, 1, 1];
        
        % set compression level
        compression_level = COMPRESSION_LEVEL;
        
    case 1
        
        % check that the matrix size is greater than 1 MB
        one_mb = 1024^2/8;
        if numel(matrix) > one_mb
            
            % set chunk size to 1 MB
            if (Nx > Ny)
                chunk_size = [one_mb, 1, 1];
            elseif (Ny > Nz)
                chunk_size = [1, one_mb, 1];
            else
                chunk_size = [1, 1, one_mb];
            end
            
            % set compression level
            compression_level = COMPRESSION_LEVEL;
            
        else
        
            % set no compression
            compression_level = 0;
            
        end
        
        
    otherwise
        error('input matrix must have 1, 2 or 3 dimensions');
    
end

% check the format of the matrix is either single precision (float in C++)
% or uint64 (unsigned long in C++)
if isa(matrix, MATRIX_DATA_TYPE_MATLAB)
    
    % set data type flags
    data_type_matlab = MATRIX_DATA_TYPE_MATLAB;
    data_type_c = MATRIX_DATA_TYPE_C;
    
elseif isa(matrix, INTEGER_DATA_TYPE_MATLAB)
    
    % set data type flags
    data_type_matlab = INTEGER_DATA_TYPE_MATLAB;
    data_type_c = INTEGER_DATA_TYPE_C;
    
else
    error('input matrix must be of type ''single'' or ''uint64''');
end

% check if the input matrix is real or complex, if complex, rearrange the
% data in the C++ format
if isreal(matrix)
    
    % set file tag
    domain_type = DOMAIN_TYPE_REAL;
    
elseif dims == 3
    
    % set file tag
    domain_type = DOMAIN_TYPE_COMPLEX;
    
    % rearrange the data so the real and imaginary parts are stored in the
    % same matrix
    matrix = cat(1,real(matrix),imag(matrix));
    matrix = reshape(matrix, Nx, 2, Ny, Nz);        
    matrix = permute(matrix, [2 1 3 4]);
    matrix = reshape(matrix, 2*Nx, Ny, Nz);
    
    % update the size of Nx
    Nx = 2*Nx;
    
elseif dims == 1
    
    % set file tag
    domain_type = DOMAIN_TYPE_COMPLEX;
    
    % rearrange the data so the real and imaginary parts are stored in the
    % same matrix
    nelems = numel(matrix);
    matrix = reshape(matrix, nelems, []);
    matrix = cat(1, real(matrix), imag(matrix));
    matrix = reshape(matrix, nelems, 2, 1, 1);        
    matrix = permute(matrix, [2 1 3 4]);
    
    % update the matrix size
    Nx = Nx * (2 - double(Nx == 1));
    Ny = Ny * (2 - double(Ny == 1));
    Nz = Nz * (2 - double(Nz == 1));
    
    % put in correct dimension
    matrix = reshape(matrix, Nx, Ny, Nz);
    
else
    
    error('currently there is no support for saving 2D complex matrices');
    
end

% allocate a holder for the new matrix within the file
if compression_level ~= 0
    
    % use compression
    h5create(filename,['/' matrix_name], [Nx, Ny, Nz], ...
                       'Datatype', data_type_matlab, ...
                       'ChunkSize', chunk_size, ...
                       'Deflate', compression_level);
                   
else
    
    % don't use compression
    h5create(filename,['/' matrix_name], [Nx, Ny, Nz], ...
                       'Datatype', data_type_matlab);
                   
end

% write the matrix into the holder within the file
h5write(filename, ['/' matrix_name], matrix);

% set attributes for the matrix (used by k-Wave++)
h5writeatt(filename, ['/' matrix_name], DOMAIN_TYPE_ATT_NAME, domain_type);
h5writeatt(filename, ['/' matrix_name], DATA_TYPE_ATT_NAME, data_type_c);