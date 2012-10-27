function saveResults( newData )
%saveResults  add some new results to the stored data
%
%   GPUBENCH.SAVERESULTS(NEWDATA) adds the supplied gpuBench result data to
%   the appropriate data-file.
%
%   Examples:
%   >> data = gpuBench;
%   >> gpubench.saveResults(data);
%
%   See also: gpuBench, gpuBenchReport

%   Copyright 2011 The MathWorks, Inc.

if ~isa( newData, 'gpubench.PerformanceData' )
    error( 'gpuBenchSaveResults:BadData', 'Input data must be a PerformanceData object as returned by gpuBench.' );
end

% Work out which data-file to use from the release
datafilename = fullfile( gpubench.getDataDir(), [newData.MATLABRelease, '.mat'] );

% Work out whether we are creating a new file or appending to an
% existing one
if exist( datafilename, 'file' ) == 2
    iAppendFile( datafilename, newData );
else
    iWriteFile( datafilename, newData );
end

end


%-------------------------------------------------------------------------%
function iAppendFile( datafilename, newResult )
% Append a new results to an existing results file. If a result for this
% card already exists it is replaced.
results = iReadFile( datafilename );
N = numel( results );

% Check to see if it already exists
deviceNames = cell(N,1);
for ii=1:N
    deviceNames{ii} = results(ii).GPUInfo.Name;
end
toRemove = strcmp( deviceNames, newResult.GPUInfo.Name );
results( toRemove ) = [];

results(end+1,1) = newResult;
iWriteFile( datafilename, results );
end

%-------------------------------------------------------------------------%
function results = iReadFile( datafilename )
% Read some existing results from file.
data = load( datafilename );
if ~isfield( data, 'results' ) || ~isa( data.results, 'gpubench.PerformanceData' )
    error( 'gpuBenchSaveResults:BadDataFile', 'Data file was corrupted, please delete it: %s', datafilename );
end
results = data.results;
end

%-------------------------------------------------------------------------%
function iWriteFile( datafilename, results ) %#ok<INUSD>
% Write the results structure to file. Currently this is just a MAT
% file, but in future we may want to go for XML to allow the use
% of the results elsewhere.
save( datafilename, 'results' );
end

