function out = gpuBenchReport( varargin )
%gpuBenchReport  create an HTML report based on some GPU performance data
%
%   gpuBenchReport(data) creates a new HTML report based on the supplied
%   PerformanceData and opens it in the system browser.
%
%   gpuBenchReport() creates an HTML report based only on the pre-stored
%   performance data.
%
%   filename = gpuBenchReport(...) returns the location of the main page
%   for the report instead of viewing it immediately.
%
%   Examples:
%   >> gpuBenchReport
%   >> data = gpuBench;
%   >> gpuBenchReport( data )
%
%   See also: gpuBench

%   Copyright 2011 The MathWorks, Inc.

narginchk( 0, 2 );
nargoutchk( 0, 1 );

% Try to get some comparison data
if nargin>0
    assert( all( cellfun( 'isclass', varargin, 'gpubench.PerformanceData' ) ) );
    [allData,allDataSummary] = getComparisonData( [varargin{:}]' );
else
    [allData,allDataSummary] = getComparisonData();
end
if isempty(allData)
    error( 'gpuBenchReport:NoData', 'No data to report' );
end

N = numel( allData );
gpubench.multiWaitbar( 'Creating GPUBench report...', 0 );

% Ensure the output folder exists
reportDir = fullfile( tempdir(), 'GPUBenchReport' );
if ~exist( reportDir, 'dir' )
    mkdir( reportDir );
end
copyFiles( reportDir );

% Create the summary page for this device
makeSummaryBarChart( reportDir, allDataSummary );
reportname = makeSummaryPage( reportDir, allDataSummary );
gpubench.multiWaitbar( 'Creating GPUBench report...', 'Increment', 1/(N+1) );

% Now create the detail pages for all devices
for ii=1:numel( allData )
    makePerformancePlots( reportDir, allData, ii );
    makeDetailPage( reportDir, allData, allDataSummary, ii );
    gpubench.multiWaitbar( 'Creating GPUBench report...', 'Increment', 1/(N+1) );
end

gpubench.multiWaitbar( 'Creating GPUBench report...', 'close' );

if nargout
    out = reportname;
else
    web( reportname );
end

end % gpuBenchReport

%-------------------------------------------------------------------------%
function [allData,allDataSummary] = getComparisonData( data )

if nargin<1 || isempty(data)
    data = [];
    % No user data, so use the current release
    thisRelease = regexp( version, 'R\d*[ab]', 'match' );
    thisRelease = thisRelease{1};
else
    % Work out which data-file to use from the release
    thisRelease = data(1).MATLABRelease;
    % Flag the user's data so that we can highlight it
    for ii=1:numel(data)
        data(ii).IsSelected = true;
    end
end

% Try to load the data for this release
otherData = loadDataForRelease( thisRelease, data );

% If no data was found, warn and try an older release
if isempty(otherData)
    allDataFiles = dir( fullfile( 'data', 'R*.mat' ) );
    allDataFiles = {allDataFiles.name};
   
    % Try each file in turn, starting from the most recent
    for ii=numel(allDataFiles):-1:1
        release = allDataFiles{ii}(1:6);
        otherData = loadDataForRelease( release, data );
        if ~isempty(otherData)
            break;
        end
    end
    % If we reach here without any data then no file had what we wanted.
    % Give up!
    if isempty( otherData )
        error( 'GPUBenchReport:NoData', ...
        ['No pre-stored data was found. Please re-download and install ', ...
        'gpuBench from the File Exchange.'] );
    else
        % Indicate that we're using an older release
        warning( 'GPUBenchReport:NoDataForRelease', ...
            ['There is no pre-stored data for this MATLAB release (%s). ' ...
            'Using data for %s instead. ', ...
            'Check the File Exchange for an updated version of gpuBench ', ...
            'that has data for this release.'], ...
            thisRelease, release );
    end
end

% Construct the summary statistics from all the results and then sort the
% original data using the summary score.
if nargin>0
    allData = [data;otherData];
else
    allData = otherData;
end

if isempty( allData )
    % No data at all, so warn and show something older
    error( 'GPUBenchReport:NoData', ...
        ['There is no pre-stored data for this MATLAB release (%s). ' ...
        'Check the File Exchange for an updated version or use gpuBench ', ...
        'to view the results for your system.'], release );
elseif isempty( otherData )
end

allDataSummary = gpubench.SummaryData( allData );
allData = allData(allDataSummary.SortOrder);
end % getComparisonData

%-------------------------------------------------------------------------%
function otherData = loadDataForRelease( release, existingData )
datafilename = fullfile( gpubench.getDataDir(), [release, '.mat'] );
if exist( datafilename, 'file' ) == 2
    otherData = load( datafilename );
    if ~isfield( otherData, 'results' ) || ~isa( otherData.results, 'gpubench.PerformanceData' )
        error( 'GPUBench:BadDataFile', 'Data file was corrupted, please delete it: %s', datafilename );
    end
    otherData = otherData.results;
    
    % We don't want to see a result for the exact same card
    if ~isempty(existingData)
        keep = true( size( otherData ) );
        for ii=1:numel( otherData )
            for jj=1:numel(existingData)
                if ~existingData(jj).IsHostData
                    keep(ii) = keep(ii) ...
                        && ~strcmpi( otherData(ii).GPUInfo.Name, existingData(jj).GPUInfo.Name );
                end
            end
            otherData(ii).IsSelected = false;
        end
        otherData = otherData(keep);
    end
else
    otherData = gpubench.PerformanceData.empty(0,1);
end
end % loadDataForRelease

%-------------------------------------------------------------------------%
function makeSummaryBarChart( outDir, summaryData )

assert( isa( summaryData, 'gpubench.SummaryData' ) );

N = numel( summaryData.DeviceName );
deviceNames = summaryData.DeviceName;
functionNames = summaryData.LongName;
flopsResults = summaryData.PeakFLOPS;
isSelected = summaryData.IsSelectedDevice;

% Sort by data-type
types = unique( summaryData.Datatype );
numOfType = zeros(size(types));
colsForType = cell(size(types));
for ii=1:numel(types)
    match = strcmp( summaryData.Datatype, types{ii} );
    numOfType(ii) = sum( match );
    colsForType{ii} = find( match );
end
colOrder = [colsForType{:}];
functionNames = functionNames(colOrder);
flopsResults = flopsResults(:,colOrder);


figh = gpubench.makeFigure( 'Summary' );
set(figh,'Position',[10 10 650 650])
bars = barh( flopsResults'/1e9 );
set( gca, 'YTickLabel', functionNames, 'YDir', 'Reverse' )
xlabel( sprintf('GFLOPS\n(higher is better)') )

% Highlight the selected entry
for ii=1:N
    if isSelected(ii)
        deviceNames{ii} = ['{\bf',deviceNames{ii},' (yours)}'];
    end
end

gpubench.legend( deviceNames{:}, 'Location', 'SouthEast' );
grid on
title( 'Performance Summary' )

% Set colors to fade from blue to red
colors = [
    linspace(0,0.6,N)
    zeros(1,N)
    linspace(0.75,0,N)
    ]';
for ii=1:N
    set( bars(ii), 'FaceColor', colors(ii,:), 'EdgeColor', 0.5*colors(ii,:) );
end

% highlight the user's result
highlightColor = [0 1 0];
selectedResults = find(summaryData.IsSelectedDevice);
N = numel(selectedResults);
for ii=1:N
    ratio = (N-ii+2)/(N+1);
    set( bars(selectedResults(ii)), ...
        'FaceColor', highlightColor*ratio, ...
        'EdgeColor', [0 0 0], ...
        'Linewidth', 2 );
end

% Add dividers between types
hold on
x = get(gca,'XLim');
for ii=1:numel(numOfType)-1
    plot( x, (numOfType(ii)+0.5)*[1 1], 'k-' )
end

gpubench.addAreaShadows( gca() );
gpubench.addGradientToAxes( gca() );

% Save the image to file for the HTML to pick up
filename = 'summarychart.png';
gpubench.captureFigure( figh, fullfile( outDir, filename ), false );
close( figh );
end % makeSummaryBarChart

%-------------------------------------------------------------------------%
function makePerformancePlots( outDir, data, thisIdx )
%Create a FLOPS plot for each results in data with the "thisIdx" result
%highlighted.

fileBase = sprintf('device%d',thisIdx);

plotNames = cell( size(data) );
for ii=1:numel(data)
    plotNames{ii} = data(ii).getDeviceName();
end
plotNames{thisIdx} = ['{\bf',plotNames{thisIdx}, ' (selected)}'];

results = data(thisIdx).Results;

for rr=1:numel( results )
    name = [results(rr).FunctionName, ' (', results(rr).DataType, ')'];
    figh = gpubench.makeFigure( name );
    color = get( gca(), 'ColorOrder' );
    
    % Plot each device's curve, if they have the relevant data
    for ii=1:numel(data)
        colorIdx = mod(ii-1,size(color,1))+1;
        % Find the corresponding result
        if hasResult( data(ii), name )
            linewidth = 1+2*(ii==thisIdx);
            gpubench.plotFLOPS( getResult( data(ii), name ), color(colorIdx,:), linewidth )
        else
            % We still need to plot something for the legend to be correct
            plot( nan, nan, 'Color', color(colorIdx,:), 'LineWidth', 1 );
        end
    end
    % Add a highlight around the peak-flops point
    colorIdx = mod(thisIdx-1,size(color,1))+1;
    thisResult = getResult( data(thisIdx), name );
    [maxVal,maxIdx] = max( 1e-9 * thisResult.NumOps ./ thisResult.Times);
    plot( thisResult.Sizes(maxIdx), maxVal, ...
        'Color', color(colorIdx,:), ...
        'Marker', 'o', ...
        'MarkerSize', 16, ...
        'Linewidth', 2 );
    
    title( name );
    gpubench.legend( plotNames{:}, 'Location', 'NorthWest' );
    gpubench.addGradientToAxes( gca() );
    
    % Due to a bug in ZBuffer renderer, the log-lines disappear when we add
    % a patch to the plot. Stick them back in now.
    xtick = str2num( get( gca, 'XTickLabel' ) ); %#ok<ST2NM>
    ylim = get( gca, 'YLim' );
    for ii=1:numel(xtick)-1
        xminortick = (2:9) * power(10,xtick(ii));
        xdata = [xminortick;xminortick;nan(1,numel(xminortick))];
        ydata = [repmat(ylim',1,numel(xminortick)); nan(1,numel(xminortick))];
        line( 'XData', xdata(:), 'YData', ydata(:), ...
            'LineStyle', ':', ...
            'Color', 0.6*[1 1 1] );
    end
    
    % Save the image to file for the HTML to pick up
    filename = [fileBase,'-',results(rr).FunctionName,'-',results(rr).DataType,'.png'];
    gpubench.captureFigure( figh, fullfile( outDir, filename ) );
    close( figh );
end
end % makePerformancePlots

%-------------------------------------------------------------------------%
function copyFiles( outDir )
% Copy the required stylesheet and other files into the report folder

dataDir = gpubench.getDataDir();

files = {
    'gpubench.css'
    'title.png'
    'background.png'
    'warning.png'
    };
for ii=1:numel(files)
    copyfile( fullfile( dataDir, files{ii} ), outDir, 'f' );
end
end % copyFiles

%-------------------------------------------------------------------------%
function reportname = makeSummaryPage( outDir, summaryData )
% Create the summary page for this device

assert( isa( summaryData, 'gpubench.SummaryData' ) );

% Find the user's data
userIdx = find( summaryData.IsSelectedDevice, 1, 'first' );
if isempty( userIdx )
    % No user data, so use the first one
    pageName = '';
else
    pageName = [': ', summaryData.DeviceName{userIdx}];
end
reportname = fullfile( outDir, 'index.html' );

fid = fopen( reportname, 'wt' );
fprintf( fid, '<html><head>\n' );
fprintf( fid, '  <title>GPU Comparison Report%s</title>\n', pageName );
fprintf( fid, '  <link rel="stylesheet" type="text/css" href="gpubench.css"/>\n' );
fprintf( fid, '  <meta author="Generated by GPUBench."/>\n' );
fprintf( fid, '</head>\n' );
fprintf( fid, '<body>\n' );

% All of the content goes in a giant table to control the width
fprintf( fid, '  <center><table class="mainlayout"><tr><td>\n' );
fprintf( fid, '  <img src="title.png"/>\n' );

% Main title
fprintf( fid, '  <h1>GPU Comparison Report%s</h1>\n', pageName );

% Summary section
fprintf( fid, '  <h2>Summary of results</h2>\n' );
% Add the description
writeBoilerPlate( fid, 'summary_intro.html' )

names = summaryData.FunctionName;
numCols = numel( names );

% Print the table header
fprintf( fid, '  <table border="0" width="100%%"><tr><td align="center">\n' );
fprintf( fid, '  <table class="summarytable" cellspacing="0">\n' );

% Titles for types
types = getDatatypes( summaryData );
colsForType = getColsForType( summaryData, types );
numOfType = cellfun( 'length', colsForType );
colOrder = [colsForType{:}];

fprintf( fid, '    <tr><th class="summarytable"></th>' );
for ii=1:numel(types)
    fprintf( fid, '      <th class="summarytable" colspan="%d">\n', ...
        numOfType(ii) );
    fprintf( fid, '        Results for data-type ''%s''<br/>(In GFLOPS)</th>\n', ...
        types{ii} );
end

% Titles for each result
fprintf( fid, '    <tr><th class="summarytable"></th>' );
for cc=1:numCols
    fprintf( fid, '<th class="summarytable">%s</th>', names{colOrder(cc)} );
end
fprintf( fid, '</tr>\n' );
% Now the body
for rr=1:numel(summaryData.DeviceName)
    fprintf( fid, '    <tr>' );
    if (summaryData.IsSelectedDevice(rr))
        nameStr = ['<b>',summaryData.DeviceName{rr},'</b>'];
        cellformat = '<td class="summarytable" align="right"><a href="device%d.html#result%u"><b>%1.2f</b></a></td>';
    else
        nameStr = summaryData.DeviceName{rr};
        cellformat = '<td class="summarytable" align="right"><a href="device%d.html#result%u">%1.2f</a></td>';
    end
    fprintf( fid, '<td class="summarytable" align="left"><a href="device%d.html">%s</a></td>', rr, nameStr );
    for cc=1:numCols
        colIdx = colOrder(cc);
        fprintf( fid, cellformat, rr, colIdx, summaryData.PeakFLOPS(rr,colIdx) / 1e9 );
    end
    fprintf( fid, '</tr>\n' );
end
fprintf( fid, '</table>\n\n' );
fprintf( fid, '<small><b>(click any device name or result to see the detailed data)</b></small><br/>\n\n' );


% Add the summary image
fprintf( fid, '<img src="summarychart.png"/>\n\n' );
fprintf( fid, '  </td></tr></table>\n' );

% Footer
fprintf( fid, '  <hr/>\n' );
writeGeneratedBy( fid );

% Close the main layout table
fprintf( fid, '  </td></tr></table></center>\n' );
fprintf( fid, '</body>\n' );
fprintf( fid, '</html>\n' );
fclose( fid );
end % makeSummaryPage

%-------------------------------------------------------------------------%
function makeDetailPage( outDir, data, summaryData, thisIdx )
% Create page of detailed information for one device
name = summaryData.DeviceName{thisIdx};
fileBase = sprintf('device%d',thisIdx);

reportname = fullfile( outDir, [fileBase,'.html'] );

fid = fopen( reportname, 'wt' );
fprintf( fid, '<html><head>\n' );
fprintf( fid, '  <title>GPU Performance Details: %s</title>\n', name );
fprintf( fid, '  <link rel="stylesheet" type="text/css" href="gpubench.css"/>\n' );
fprintf( fid, '  <meta author="Generated by GPUBench."/>\n' );
fprintf( fid, '</head>\n' );
fprintf( fid, '<body>\n' );

% All of the content goes in a giant table to control the width
fprintf( fid, '  <center><table class="mainlayout"><tr><td>\n' );
fprintf( fid, '  <a href="index.html"><img src="title.png" border="0"/></a>\n' );

% Main title
fprintf( fid, '  <h1>GPU Performance Details: %s</h1>\n', name );

% Contents section
fprintf( fid, '  <table cellspacing="10"><tr><td valign="top">\n' );
fprintf( fid, '    <b>Contents:</b>\n' );
fprintf( fid, '  </td><td valign="top">\n' );
names = summaryData.LongName;
fprintf( fid, '    <ul>\n' );
fprintf( fid, '      <li><a href="#config">System Configuration</a></li>\n' );

% Sort the contents by type
types = getDatatypes( summaryData );

for tt=1:numel( types )
    fprintf( fid, '      <li>Results for datatype %s</a><ul>\n', types{tt} );
    colsForType = getColsForType( summaryData, types{tt} );
    for nn=1:numel( colsForType )
        myCol = colsForType(nn);
        fprintf( fid, '        <li><a href="#result%u">%s</a></li>\n', myCol, names{myCol} );
    end
    fprintf( fid, '      </ul></li>\n' );
end
fprintf( fid, '    </ul>\n' );
fprintf( fid, '  </tr></table>\n' );
fprintf( fid, '  <br/>\n' );


% Add a section showing the operating environment
fprintf( fid, '  <a name="config"><h2>System Configuration</h2></a>\n' );

% If this isn't the user's system, make that clear
if ~data(thisIdx).IsSelected
    fprintf( fid, '  <table border="0"><tr>\n' );
    fprintf( fid, '    <td valign="middle"><img src="warning.png"/></td>\n' );
    fprintf( fid, '    <td valign="baseline"><b><font color="#660000">Note that this is previously stored data and does not reflect your system configuration.</font></b></td>\n' );
    fprintf( fid, '  </tr></table>\n' );
end

fprintf( fid, '  <p><b>MATLAB Release:</b> %s</p>\n', data(thisIdx).MATLABRelease );
fprintf( fid, '  <table cellspacing="10"><tr><td valign="top">\n' );
fprintf( fid, '    <p><b>Host</b></p>\n' );
writeStructTable( fid, data(thisIdx).CPUInfo );
% - only add the GPU if we ran on it
if ~data(thisIdx).IsHostData
    fprintf( fid, '  </td><td valign="top">\n' );
    fprintf( fid, '    <p><b>GPU</b></p>\n' );
    writeStructTable( fid, data(thisIdx).GPUInfo );
end
fprintf( fid, '  </td></tr></table>\n' );
fprintf( fid, '  <br/>\n' );

% Add one section per result
names = summaryData.LongName;
for nn=1:numel( names )
    fprintf( fid, '  <a name="result%u"><h2>Results for %s</h2></a>\n', nn, names{nn} );
    if ~hasResult( data(thisIdx), names{nn} )
        % No results for this function
        fprintf( fid, '  <p>No results found for %s.</p>\n', names{nn} );
        continue;
    end
    myResult = getResult( data(thisIdx), names{nn} );
    % See if there's a description for this function
    writeBoilerPlate( fid, [myResult.FunctionName,'.html'] )
    
    fprintf( fid, '  <table cellspacing="0"><tr><td valign="top">\n' );
    
    fprintf( fid, '    <table><tr><td>\n');
    fprintf( fid, '      <b>Raw data for %s - %s</b>\n', name, names{nn} );
    fprintf( fid, '    </td></tr><tr><td>\n');
    
    % Print the table header
    fprintf( fid, '    <table class="summarytable" cellspacing="0">\n' );
    fprintf( fid, '      <tr>' );
    fprintf( fid, '<th class="summarytable">Array size<br/>(elements)</th>' );
    fprintf( fid, '<th class="summarytable">Num<br/>Operations</th>' );
    fprintf( fid, '<th class="summarytable">Time<br/>(ms)</th>' );
    fprintf( fid, '<th class="summarytable">GigaFLOPS</th>' );
    fprintf( fid, '</tr>\n' );
    % Now one row per size
    sizes = myResult.Sizes;
    flops = myResult.NumOps;
    times = myResult.Times;
    [~,peakIdx] = max( flops ./ times );
    baseFormatStr1 = '<td class="summarytable" align="right">';
    baseFormatStr2 = '</td>';
    for ss=1:numel(sizes)
        % Highlight the peak FLOPS row
        if ss==peakIdx
            formatStr1 = [baseFormatStr1,'<font color="#0000dd">'];
            formatStr2 = ['</font>',baseFormatStr2];
        else
            formatStr1 = baseFormatStr1;
            formatStr2 = baseFormatStr2;
        end
        fprintf( fid, '      <tr>' );
        fprintf( fid, [formatStr1,'%s',formatStr2], num2strWithCommas(sizes(ss)) );
        fprintf( fid, [formatStr1,'%s',formatStr2], num2strWithCommas(flops(ss)) );
        fprintf( fid, [formatStr1,'%2.2f',formatStr2], times(ss)*1000 );
        fprintf( fid, [formatStr1,'%2.2f',formatStr2], flops(ss)/times(ss)/1e9 );
        fprintf( fid, '</tr>\n' );
    end
    fprintf( fid, '    </table>\n' );
    fprintf( fid, '<center><small>(<code>N</code> gigaflops = <code>Nx10<sup>9</sup></code> operations per second)</small></center><br/>\n\n' );
    
    fprintf( fid, '    </td></tr></table>\n');
    
    % Add the image
    fprintf( fid, '  </td><td valign="top">\n' );
    fprintf( fid, '    <img src="%s-%s-%s.png"/>\n', ...
        fileBase, myResult.FunctionName, myResult.DataType );
    fprintf( fid, '  </td></tr></table>\n' );
end

% Footer
fprintf( fid, '  <hr/>\n' );

fprintf( fid, '  <table width="100%%"><tr><td align="left">\n' );
writeGeneratedBy( fid );
fprintf( fid, '  </td><td align="right">\n' );
fprintf( fid, '    <small><a href="index.html"><i>Back to summary</i></a></small>\n' );
fprintf( fid, '  </td></tr></table>\n' );

% Close the main layout table
fprintf( fid, '  </td></tr></table></center>\n' );
fprintf( fid, '</body>\n' );
fprintf( fid, '</html>\n' );
fclose( fid );
end % makeDetailPage

%-------------------------------------------------------------------------%
function writeStructTable( fid, data )
assert( isstruct( data ) && isscalar( data ) );
fprintf( fid, '    <table class="summarytable" cellspacing="0">\n' );
fields = fieldnames( data );
for ff=1:numel( fields )
    fprintf( fid, '      <tr><th class="summarytable" align="left">%s</th>', fields{ff} );
    fprintf( fid, '<td class="summarytable">' );
    x = data.(fields{ff});
    if ischar( x )
        fprintf( fid, '%s', x );
    elseif isinteger( x )
        fprintf( fid, '%d', x );
    else
        % Try to let MATLAB do the right thing
        fprintf( fid, '%g', x );
    end
    fprintf( fid, '</td></tr>\n' );
end
fprintf( fid, '    </table>\n' );
end % writeStructTable

%-------------------------------------------------------------------------%
function writeBoilerPlate( outFid, filename )
%Read some boiler-plate HTML and paste it into the supplied output file
filename = fullfile( gpubench.getDataDir(), filename );
inFid = fopen( filename, 'rt' );
if inFid<=0
    warning( 'gpuBenchReport:MissingBoilerPlateFile', ...
        'Input file could not be opened: %s', filename );
    return;
end
txt = fread( inFid );
fwrite( outFid, txt );
fclose( inFid );
end % writeBoilerPlate

%-------------------------------------------------------------------------%
function writeGeneratedBy( outFid )
%Write the "generated by" string into the footer

fprintf( outFid, '<small><i>Generated by gpuBench v%s: %s</i></small>\n', ...
    gpubench.version(), datestr( now(), 'yyyy-mm-dd HH:MM:SS' ) );
end % writeGeneratedBy

%-------------------------------------------------------------------------%
function str = num2strWithCommas( num )
%Convert an integer into a string with commas separating sets of 3 digits
%
%  e.g. num2StrWithCommas(12345678) = '12,345,678'

% First convert using the standard method
baseStr = num2str( abs(num) );
% now insert some commas.
% pad to a multiple of 3
padding = 3 - (mod(length(baseStr)-1,3)+1);
str = [repmat(' ',1,padding), baseStr];
numCols = length(str)/3;
str = [reshape(str,3,numCols);repmat(',',1,numCols)];
str = strtrim( str(1:end-1) );
% Finally, re-insert the sign
if num<0
    str = ['-',str];
end
end % num2StrWithCommas