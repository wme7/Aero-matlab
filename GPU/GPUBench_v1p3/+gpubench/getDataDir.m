function dataDir = getDataDir()
%getDataDir  get the full path to the folder containing the GPUBench data

%   Copyright 2011 The MathWorks, Inc.

dataDir = fullfile( fileparts( fileparts( mfilename( 'fullpath' ) ) ), 'data' );

