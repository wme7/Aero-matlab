function [outGPU,outHost] = gpuBench()
%GPUBENCH  MATLAB GPU Benchmark
%   GPUBENCH times different MATLAB GPU tasks and compares the execution
%   speed with the speed of several other GPUs.  The tasks are:
%
%    Backslash   Matrix left-division.    Floating point, regular memory access.
%    MTimes      Matrix multiplication.   Floating point, regular memory access.
%    FFT         Fast Fourier Transform.  Floating point, irregular memory access.
%
%   Each task is run for a range of array sizes and the results are tabulated
%   in an HTML report.  GPUBENCH can take several minutes to complete - please
%   be patient! Note that if your GPU is also driving your monitor then
%   the display may become unresponsive during testing.
%
%   GPUBENCH runs each of the tasks and shows a report indicating how the
%   current GPU compares to other systems.
%
%   T = GPUBENCH returns a data structure containing all of the results and
%   does not generate the report.
%
%   Fluctuations of up to ten percent in the measured times of repeated
%   runs on a single machine are not uncommon.  Your own mileage may vary.
%
%   This benchmark is intended to compare performance different GPUs on one
%   particular version of MATLAB.  It does not offer direct comparisons
%   between different versions of MATLAB.
%
%   See also: BENCH, gpuBenchReport

% Unused tasks:
%    Mandelbrot  Calculate a Mandelbrot Set.  Floating point, regular memory access.

%   Author: Ben Tordoff
%   Copyright 2011-2012 The MathWorks, Inc.

% Check for the right MATLAB version and availability of PCT
gpubench.checkMATLABVersion();
gpubench.checkPCT();

% Check for a GPU. We give the option of running without a GPU so that
% users can evaluate what benefits a GPU might give.
hasGPU = parallel.gpu.GPUDevice.isAvailable();
if ~hasGPU
    title = 'Continue without a GPU?';
    question = ['The GPU could not be used. ' ...
        'Do you wish to continue and collect results for your CPU?'];
    buttons = {'Collect CPU results', 'Stop'};
    answer = questdlg(question, title, buttons{:}, buttons{end});
    if ~strcmp(answer,buttons{1})
        warning( 'GPUBench:NoGPU', 'No GPU was available for GPUBench to use.' );
        return;
    end
end

% Initialize the data object
release = regexp( version, 'R\d*[ab]', 'match' );
gpuData = gpubench.PerformanceData( ...
    release{1}, ...
    gpubench.cpuinfo(), ...
    gpubench.gpuinfo(), ...
    now() );
hostData = gpuData;
hostData.IsHostData = true;

% Do we need to measure the host stuff?
doHost = (nargout~=1);
numTasks = 6*(hasGPU+doHost);
reps = 3;
progressTitle = 'Running GPUBench...';
gpubench.multiWaitbar( progressTitle, 0 );

if hasGPU
    gpuData = runBackslash( gpuData, reps, 'single', 'GPU', progressTitle, numTasks );
    gpuData = runBackslash( gpuData, reps, 'double', 'GPU', progressTitle, numTasks );

    gpuData = runMTimes( gpuData, reps, 'single', 'GPU', progressTitle, numTasks );
    gpuData = runMTimes( gpuData, reps, 'double', 'GPU', progressTitle, numTasks );

    gpuData = runFFT( gpuData, reps, 'single', 'GPU', progressTitle, numTasks );
    gpuData = runFFT( gpuData, reps, 'double', 'GPU', progressTitle, numTasks );

    % gpuData = runMandelbrot( gpuData, reps, 'double', 'GPU', progressTitle, numTasks );
end

if doHost
    hostData = runBackslash( hostData, reps, 'single', 'Host', progressTitle, numTasks );
    hostData = runBackslash( hostData, reps, 'double', 'Host', progressTitle, numTasks );
    
    hostData = runMTimes( hostData, reps, 'single', 'Host', progressTitle, numTasks );
    hostData = runMTimes( hostData, reps, 'double', 'Host', progressTitle, numTasks );
    
    hostData = runFFT( hostData, reps, 'single', 'Host', progressTitle, numTasks );
    hostData = runFFT( hostData, reps, 'double', 'Host', progressTitle, numTasks );
    
    % hostData = runMandelbrot( hostData, reps, 'double', 'Host', progressTitle, numTasks );
end

gpubench.multiWaitbar( progressTitle, 'Close' );

if nargout
    % User requested raw data
    outGPU = gpuData;
    outHost = hostData;
else
    % Produce report
    reportData = {};
    if hasGPU
        reportData{end+1} = gpuData;
    end
    if doHost
        reportData{end+1} = hostData;
    end
    web( gpuBenchReport( reportData{:} ) );
end


%-------------------------------------------------------------------------%
function data = runFFT( data, reps, type, device, mainProgressTitle, numTasks )
% Work out the maximum size we should run
safetyFactor = 6;
sizes = getTestSizes( type, safetyFactor, device );
times = inf( size( sizes ) );
worstTime = 0;

progressTitle = sprintf( 'FFT (%s, %s)', device, type );
progressTotal = sum(sizes);
gpubench.multiWaitbar( progressTitle, 0 );

for ii=1:numel(sizes)
    % Check for getting close to time-out
    if tooCloseToTimeout( worstTime, device )
        fprintf( 'Skipping FFT of size %u to prevent timeout.\n', sizes(ii) );
        times(ii) = nan;
        continue;
    end        
    N = sizes(ii);
    try
        A = complex( rand( N, 1, type ), rand( N, 1, type ) );
        if strcmpi( device, 'GPU' )
            A = gpuArray(A);
        end
        
        for rr=1:reps
            t = tic();
            B = fft(A); %#ok<NASGU>
            elapsedTime = gtoc(t);
            times(ii) = min( times(ii), elapsedTime );
            worstTime = max( worstTime, elapsedTime );
            clear B;
            % Update both progress bars
            inc = sizes(ii)/(reps*progressTotal);
            gpubench.multiWaitbar( progressTitle, 'Increment', inc );
            gpubench.multiWaitbar( mainProgressTitle, 'Increment', inc/numTasks );
        end
    catch err %#ok<NASGU>
        fprintf( 'discarded FFT of size %u.\n', N );
        times(ii) = nan;
    end
end
gpubench.multiWaitbar( progressTitle, 'Close' );

% Clear any dud results
sizes(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'FFT', type, sizes, 5*sizes.*log2(sizes), times );


%-------------------------------------------------------------------------%
function data = runMTimes( data, reps, type, device, mainProgressTitle, numTasks )
safetyFactor = 2;
sizes = getTestSizes( type, safetyFactor, device );

times = inf( size( sizes ) );
worstTime = 0;

progressTitle = sprintf( 'MTimes (%s, %s)', device, type );
progressTotal = sum(sizes);
gpubench.multiWaitbar( progressTitle, 0 );

N = round( sqrt( sizes ) );
for ii=1:numel(sizes)
    % Check for getting close to time-out
    if tooCloseToTimeout( worstTime, device )
        fprintf( 'Skipping MTimes of %ux%u to prevent timeout.\n', N(ii), N(ii) );
        times(ii) = nan;
        continue;
    end        
    
    try
        A = rand( N(ii), N(ii), type );
        B = rand( N(ii), N(ii), type );
        if strcmpi( device, 'GPU' )
            A = gpuArray(A);
            B = gpuArray(B);
        end
        for rr=1:reps
            t = tic();
            C = A*B; %#ok<NASGU>
            elapsedTime = gtoc(t);
            times(ii) = min( times(ii), elapsedTime );
            worstTime = max( worstTime, elapsedTime );
            clear C;
            % Update both progress bars
            inc = sizes(ii)/(reps*progressTotal);
            gpubench.multiWaitbar( progressTitle, 'Increment', inc );
            gpubench.multiWaitbar( mainProgressTitle, 'Increment', inc/numTasks );
        end
    catch err %#ok<NASGU>
        fprintf( 'discarded MTimes of %ux%u.\n', N(ii), N(ii) );
        times(ii) = nan;
    end
end
gpubench.multiWaitbar( progressTitle, 'Close' );

% Clear any dud results
N(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'MTimes', type, N.*N, N.*N.*(2.*N-1), times );



%-------------------------------------------------------------------------%
function data = runBackslash( data, reps, type, device, mainProgressTitle, numTasks )
safetyFactor = 1.5;
sizes = getTestSizes( type, safetyFactor, device );

% Limit the sizes to 1e8 for now to prevent problems
sizes(sizes>1e8) = [];

times = inf( size( sizes ) );
worstTime = 0;

progressTitle = sprintf( 'Backslash (%s, %s)', device, type );
progressTotal = sum(sizes);
gpubench.multiWaitbar( progressTitle, 0 );

N = round( sqrt( sizes ) );
for ii=1:numel(sizes)
    % Check for getting close to time-out
    if tooCloseToTimeout( worstTime, device )
        fprintf( 'Skipping Backslash of %ux%u to prevent timeout.\n', N(ii), N(ii) );
        times(ii) = nan;
        continue;
    end        
    try
        A = 10*eye( N(ii), N(ii), type ) + rand( N(ii), N(ii), type );
        b = rand( N(ii), 1, type );
        if strcmpi( device, 'GPU' )
            A = gpuArray(A);
            b = gpuArray(b);
        end
        for rr=1:reps
            t = tic();
            C = A\b; %#ok<NASGU>
            elapsedTime = gtoc(t);
            times(ii) = min( times(ii), elapsedTime );
            worstTime = max( worstTime, elapsedTime );
            clear C;
            % Update both progress bars
            inc = sizes(ii)/(reps*progressTotal);
            gpubench.multiWaitbar( progressTitle, 'Increment', inc );
            gpubench.multiWaitbar( mainProgressTitle, 'Increment', inc/numTasks );
        end

    catch err %#ok<NASGU>
        fprintf( 'discarded Backslash of %ux%u.\n', N(ii), N(ii) );
        times(ii) = nan;
    end
end
gpubench.multiWaitbar( progressTitle, 'Close' );

% Clear any dud results
N(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'Backslash', type, N.*N, round(2/3*N.^3 + 3/2*N.^2), times );


%-------------------------------------------------------------------------%
function data = runMandelbrot( data, reps, type, device, mainProgressTitle, numTasks ) %#ok<DEFNU>
% This task will only run on the GPU
safetyFactor = 3;
sizes = getTestSizes( type, safetyFactor, device );

times = inf( size( sizes ) );
worstTime = 0;
numops = inf( size( sizes ) );
maxIterations = 200;
xlim = [-2, 0.5];
ylim = [ -1.25,  1.25];

progressTitle = sprintf( 'Mandelbrot (%s, %s)', type, device );
progressTotal = sum(sizes);
gpubench.multiWaitbar( progressTitle, 0 );

for ii=1:numel(sizes)
    gridSize = round( sqrt( sizes(ii) ) );
    % Check for getting close to time-out
    if tooCloseToTimeout( worstTime, device )
        fprintf( 'Skipping Mandelbrot of size %ux%u to prevent timeout.\n', gridSize, gridSize );
        times(ii) = nan;
        continue;
    end        
    if strcmpi( device, 'GPU' )
        try
            
            x = parallel.gpu.GPUArray.linspace( xlim(1), xlim(2), gridSize );
            y = parallel.gpu.GPUArray.linspace( ylim(1), ylim(2), gridSize );
            [xGrid,yGrid] = meshgrid( x, y );
            
            % Calculate
            for rr=1:reps
                t = tic();
                count = arrayfun( @processMandelbrotElement, xGrid, yGrid, maxIterations );
                elapsedTime = gtoc(t);
                times(ii) = min( times(ii), elapsedTime );
                worstTime = max( worstTime, elapsedTime );
            end
            % Use the count to work out the number of operations
            % Each iteration of a single element requires:
            %  * abs(complex) = 3 flop
            %  * count+1 = 1 flop
            %  * z*z + z0 = 8 flop
            numops(ii) = gather( sum(count(:)*12) );
            
            clear count;
            % Update both progress bars
            inc = sizes(ii)/(reps*progressTotal);
            gpubench.multiWaitbar( progressTitle, 'Increment', inc );
            gpubench.multiWaitbar( mainProgressTitle, 'Increment', inc/numTasks );
            
        catch err %#ok<NASGU>
            fprintf( 'discarded Mandelbrot of %ux%u.\n', gridSize, gridSize );
            times(ii) = nan;
        end
        
    else
        % Host version. This takes too long to do several repeats so we
        % just run it once.
        x = linspace( xlim(1), xlim(2), gridSize );
        y = linspace( ylim(1), ylim(2), gridSize );
        [xGrid,yGrid] = meshgrid( x, y );
        z0 = complex(xGrid,yGrid);
        t = tic();
        z = z0;
        count = zeros(size(z));
        for n = 1:maxIterations
            inside = ((real(z).^2 + imag(z).^2) <= 4);
            count = count + inside;
            z = z.*z + z0;
        end
        times(ii) = toc(t);
        % Each iteration of a single element requires:
        %  * inside check = 3 flop
        %  * count+1 = 1 flop
        %  * z*z + z0 = 8 flop
        % Since every element does the same amount of work in this
        % version, the operation count is simply 12*numel*maxIters
        numops(ii) = 12*numel(z)*maxIterations;
        
        % Update both progress bars
        gpubench.multiWaitbar( progressTitle, 'Increment', sizes(ii)/progressTotal );
        gpubench.multiWaitbar( mainProgressTitle, 'Increment', (sizes(ii)/progressTotal)/numTasks );
    end
end
gpubench.multiWaitbar( progressTitle, 'Close' );

% Clear any dud results
sizes(isnan( times )) = [];
numops(isnan( times )) = [];
times(isnan( times )) = [];

data = addResult( data, 'Mandelbrot', type, sizes, numops, times );

%-------------------------------------------------------------------------%
function elapsedTime = gtoc( timer )
% Wait for GPU operations to complete the call toc
persistent hasWait;
if isempty(hasWait)
    try
        wait(gpuDevice);
        hasWait = true;
    catch err %#ok<NASGU>
        hasWait = false;
    end
elseif hasWait
    wait(gpuDevice);
end
elapsedTime = toc(timer);

%-------------------------------------------------------------------------%
function sizes = getTestSizes( type, safetyFactor, device )
% Return the maximum number of elements that will fit in the device memory
elementSize = gpubench.sizeof( type );
if strcmpi( device, 'Host' )
    % On the host everything takes longer, so don't go as far
    safetyFactor = safetyFactor*2;
end

% Use as much memory as we can.
if parallel.gpu.GPUDevice.isAvailable()
    gpu = gpuDevice();
    freeMem = gpu.FreeMemory;
else
    % No GPU to get memory size, so just go for 4GB
    freeMem = 4*2^30;
end
maxNumElements = floor( freeMem / (elementSize*safetyFactor) );
if isnan( maxNumElements ) || maxNumElements < 1e6
    error( 'gpuBench:NotEnoughMemory', 'Not enough free device memory to run tasks' );
end

% We want powers of two up to this size
maxPower = floor( log2( maxNumElements ) );
sizes = power( 2, 10:2:maxPower );


%-------------------------------------------------------------------------%
function stopNow = tooCloseToTimeout( time, device )
% Shoulkd a test stop early to avoid triggering the device time-out?
stopNow = false;
if strcmpi( device, 'Host' )
    % On the host there is no time limit
else
    gpu = gpuDevice();
    % If the kernel has a timeout it is typically 2-5 seconds. If we have
    % just done a size that takes more than a quarter of a second, the next
    % size will likely trigger the timeout.
    stopNow = (gpu.KernelExecutionTimeout && time>0.25);
end



%-------------------------------------------------------------------------%
function count = processMandelbrotElement(x0,y0,maxIterations)
% Evaluate the Mandelbrot function for a single element
%
%   m = processMandelbrotElement(x0,y0,maxIterations) evaluates the
%   number of steps before the complex value (x0,y0) jumps outside a circle
%   of radius two on the complex plane. Each iteration involves mapping
%   z=z^2+z0 where z0=x0+i*y0. The return value is the log of the
%   iteration count at escape or maxIterations if the point did not escape.
z0 = complex(x0,y0);
z = z0;
count = 0;
while (count < maxIterations) ...
        && ((real(z)*real(z) + imag(z)*imag(z)) <= 4)
    count = count + 1;
    z = z*z + z0;
end
