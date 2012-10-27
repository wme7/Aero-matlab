function window = mandelbrotViewer()
%mandelbrotViewer  view and explore the Mandelbrot set using a GPU
%
%   mandelbrotViewer() opens a MATLAB figure window showing the Mandelbrot
%   set. Use the usual zoom and pan controls from the figure window toolbar
%   to navigate around and explore, or click "animate" to see a pre-defined
%   path through the set. You can move back to the initial view at any time
%   by clicking the "reset" button or add the current view to the animation
%   list using "add".
%
%   The control panel can be hidden using the right-hand toolbar button.
%
%   A selector allows a choice of four ways to calculate each frame:
%
%   1. CPU: All calculations are performed by MATLAB on the host CPU. The
%   algorithm is fully vectorized and avoids indexing to give an efficient
%   calculation. Even so, this may take a few seconds to calculate each
%   frame.
%
%   2. GPU (simple): The same algorithm as (1) is used but the input
%   coordinates are switched to being stored on the GPU. This causes
%   MATLAB to operate on the resulting data array on the GPU with no other
%   code change. This gives some speedup at virtually no coding cost.
%
%   3. GPU Arrayfun: Now we change the code so that instead of many
%   operations running on the full data matrix, we now specify a single
%   "calculateElement" operation and run it for each element. Roughly
%   speaking, MATLAB translates one element into one thread on the GPU,
%   giving huge speedups for large arrays. Note that all code is still
%   written in MATLAB and the user needs no knowledge of how GPU kernels
%   are constructed and executed.
%
%   4. CUDAKernel: Taking things to the limit, we now hand-craft the
%   element-wise algorithm used in (3) in CUDA C++. The resulting kernel is
%   called from MATLAB using the CUDAKernel system, requiring the user to
%   specify the thread and block arrangements to use.
%
%   Note that version 3 gets us most of the speedup achieved by the
%   hand-crafted CUDA (version 4) but without any need to leave the comfort
%   of MATLAB!
%
%   See also:  gpuArray, mandelbrotViewerProcessElement

%   Copyright 2010-2011 The Mathworks, Inc.

% Check that we are running in R2011a or above and have a GPU
matlabVersionCheck();
gpuCheck();

% Define some global (to this file) data structures so that they can be
% used by all the helper functions.
data = createData();
gui = createGUI();
% Make sure the image is updated now that the window is onscreen
redraw();

% Return the window handle if requested
if nargout
    window = gui.Window;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here are the four different Mandelbrot Set computations

    function logCount = computeMandelbrotCPU( xlim, numx, ylim, numy, maxIters )
        % Create a view of the Mandelbrot set using only the CPU.
        % This is the base-line version of the algorithm that is adapted to
        % run on the GPU in different ways in the functions below.
        
        % Create the input arrays
        escapeRadius2 = 400; % Square of escape radius
        x = linspace( xlim(1),  xlim(2), numx );
        y = linspace( ylim(1),  ylim(2), numy );
        [x0,y0] = meshgrid(x, y);
        count = zeros( size( x0 ) );
        z0 = complex( x0, y0 );
        z = z0;
        
        % Calculate
        for n = 0:maxIters
            inside = ((real(z).^2 + imag(z).^2) <= escapeRadius2);
            count = count + inside;
            z = inside.*(z.*z + z0) + (1-inside).*z;
        end
        magZ2 = real(z).^2 + imag(z).^2;
        logCount = log( count + 1 - log( log( max(magZ2,escapeRadius2) ) / 2 ) / log(2) );
    end % computeMandelbrotCPU


    function logCount = computeMandelbrotGPU( xlim, numx, ylim, numy, maxIters )
        % Compute using GPUArray overloads.
        % In this version the main calculation is exactly as it was for the
        % CPU version, we have simply changed the input grid to be on the
        % GPU. When MATLAB encounters GPU data it tries to run any
        % functions on the GPU. This provides a simple way to see if the
        % GPU helps without altering your code.
        
        % Setup the input grid on the GPU
        escapeRadius2 = 400; % Square of escape radius
        x = parallel.gpu.GPUArray.linspace( xlim(1),  xlim(2), numx );
        y = parallel.gpu.GPUArray.linspace( ylim(1),  ylim(2), numy );
        count = parallel.gpu.GPUArray.zeros( numy, numx );
        [x0,y0] = meshgrid(x, y);
        z0 = complex( x0, y0 );
        
        % Calculate
        z = z0;
        for n = 0:maxIters
            inside = ((real(z).^2 + imag(z).^2) <= escapeRadius2);
            count = count + inside;
            z = inside.*(z.*z + z0) + (1-inside).*z;
        end
        magZ2 = real(z).^2 + imag(z).^2;
        logCount = log( count + 1 - log( log( max(magZ2,escapeRadius2) ) / 2 ) / log(2) );
        
        % Gather the result back to the CPU
        logCount = gather( logCount );
    end % computeMandelbrotGPU

    function logCount = computeMandelbrotArrayFun( xlim, numx, ylim, numy, maxIters )
        % Compute using GPU arrayfun.
        % The second way in which MATLAB can use the GPU is by placing your
        % algorithm inside a helper function and calling it using ARRAYFUN
        % with some GPU data as input. The helper function is converted
        % into native GPU code (PTX) and each element of the input array is
        % processed in a separate GPU thread. The helper function must
        % operate only one scalars. Here we have taken the "calculate" code
        % above and put it in a helper "mandelbrotViewerProcessElement"
        % that will be converted into native code and run on the GPU.
        
        % Create the input arrays
        escapeRadius2 = 400;
        x = parallel.gpu.GPUArray.linspace( xlim(1),  xlim(2), numx );
        y = parallel.gpu.GPUArray.linspace( ylim(1),  ylim(2), numy );
        [x0,y0] = meshgrid(x, y);
        
        % Calculate
        [logCount] = arrayfun( @mandelbrotViewerProcessElement, x0, y0, ...
            escapeRadius2, maxIters );
        
        % Gather the result back to the CPU
        logCount = gather( logCount );
    end % computeMandelbrotCPU

    function logCount = computeMandelbrotCUDAKernel( xlim, numx, ylim, numy, maxIters )
        % Use pre-existing CUDA/C++ code.
        % The final way in which MATLAB can use the GPU is by calling some
        % hand-written CUDA code. The "CUDAKernel" interface allows the
        % function to be specified along with the number of threads and
        % blocks to use. This requires some knowledge of how GPUs work, but
        % does allow you to easily use existing CUDA kernels with MATLAB
        % data.
        
        % Create the input arrays
        escapeRadius = 20;
        x = parallel.gpu.GPUArray.linspace( xlim(1),  xlim(2), numx );
        y = parallel.gpu.GPUArray.linspace( ylim(1),  ylim(2), numy );
        [x0,y0] = meshgrid(x, y);
        
        % Make sure we have sufficient blocks to cover the whole array
        numElements = numel( x0 );
        data.Kernel.ThreadBlockSize = [data.Kernel.MaxThreadsPerBlock,1,1];
        data.Kernel.GridSize = [ceil(numElements/data.Kernel.MaxThreadsPerBlock),1];
        
        % Call the kernel
        logCount = parallel.gpu.GPUArray.zeros( size( x0 ) );
        logCount = feval( data.Kernel, logCount, ...
            x0, y0, ...
            escapeRadius, maxIters, numElements );
        logCount = gather( logCount );
    end % computeMandelbrotCUDAKernel




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Everything else is a callback or helper function.

    function out = createData()
        out = struct( ...
            'MaxIterations', 5000, ...
            'OrigXLim', [-2 1], ...
            'OrigY', 0, ...
            'XLim', [-2 1], ...
            'Y', 0, ...
            'CalculationMethods', {{
            'CPU'
            'GPU (simple)'
            'GPU ArrayFun'
            'CUDAKernel'
            }}, ...
            'IsAnimating', false, ...
            'NextLocation', 1, ...
            'LastFrameTime', now(), ...
            'WindowPixelSize', [100 100], ...
            'SelectedCalculationMethod', 'CUDAKernel', ...
            'ControlsVisible', true, ...
            'WriteVideo', false, ...
            'VideoWriter',  [] );
        
        % open writer
        if out.WriteVideo
            out.VideoWriter = VideoWriter('out.avi');
            out.VideoWriter.FrameRate = 20;
            out.VideoWriter.Quality = 90;
            open( out.VideoWriter );
        end
        
        % Read the location list from a local file
        out.LocationList = readLocationList();
        
        % Load the CUDA kernel
        out.Kernel = loadKernel();
    end % createData

    function out = createGUI()
        % Create the GUI, storing handles in the global GUI structure
        out.Window = figure( ...
            'Name', 'Mandelbrot viewer v1.1', ...
            'NumberTitle', 'off', ...
            'HandleVisibility', 'off', ...
            'MenuBar', 'none', ...
            'ToolBar', 'figure', ...
            'Renderer', 'ZBuffer' ); % Can't use painters as colormaps are broken for >256 colors!
        out.MainAxes = axes( ...
            'Parent', out.Window, ...
            'Position', [0 0 1 1], ...
            'XLim', data.XLim, ...
            'YLim', [-1 1], ...
            'CLim', log([1 data.MaxIterations]), ...
            'XTick', [], 'YTick', [], ...
            'DataAspectRatio', [1 1 1] );
        out.Image = image( ...
            'XData', [0 1], ...
            'YData', [0 1], ...
            'XLimInclude', 'off', ...
            'YLimInclude', 'off', ...
            'CData', nan, ...
            'CDataMapping', 'Scaled', ...
            'HandleVisibility', 'off', ...
            'Parent', out.MainAxes );
        % Add a line so that zooming works. Strange but true.
        line( 'Parent', out.MainAxes, 'XData', [-2 2], 'YData', [-2 2], ...
            'Visible', 'off', ...
            'HitTest', 'off' );
        set( out.MainAxes, 'XLimMode', 'manual', 'YLimMode', 'manual', ...
            'CLim', [0 1] );
        colormap( out.MainAxes, jet2(1000) );
        
        out.ControlPanel = uipanel( ...
            'Parent', out.Window, ...
            'BackgroundColor', 'k', ...
            'Units', 'Pixels', ...
            'Position', [10 10 155 78] );
        
        % Some text for showing compute time
        out.ComputeText = uicontrol( ...
            'Style', 'Text', ...
            'String', 'Computed in 0ms', ...
            'BackgroundColor', 'k', ...
            'ForegroundColor', 'g', ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'Position', [5 16 145 14] );
        out.FrameRateText = uicontrol( ...
            'Style', 'Text', ...
            'String', 'Displaying at 0fps', ...
            'BackgroundColor', 'k', ...
            'ForegroundColor', 'g', ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'Position', [5 2 145 14] );
        
        % Create a drop-down for selecting the calculation method
        out.MethodSelector = uicontrol( ...
            'Style', 'PopupMenu', ...
            'String', data.CalculationMethods, ...
            'Value', numel( data.CalculationMethods ), ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'BackgroundColor', 0.8*[1 1 1], ...
            'Position', [5 54 145 16], ...
            'Callback', @onCalculationMethodChanged );
        out.AddButton = uicontrol( ...
            'Style', 'ToggleButton', ...
            'String', 'Add', ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'BackgroundColor', 0.6*[1 1 1], ...
            'Position', [5 30 45 16], ...
            'TooltipString', 'Add the current location to the animation list', ...
            'Callback', @onAddPressed );
        out.ResetButton = uicontrol( ...
            'Style', 'ToggleButton', ...
            'String', 'Reset', ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'BackgroundColor', 0.6*[1 1 1], ...
            'Position', [55 30 45 16], ...
            'TooltipString', 'Reset the view to the top', ...
            'Callback', @onResetPressed );
        out.AnimateButton = uicontrol( ...
            'Style', 'ToggleButton', ...
            'String', 'Animate', ...
            'FontSize', 7, ...
            'Parent', out.ControlPanel, ...
            'BackgroundColor', 0.6*[1 1 1], ...
            'Position', [105 30 45 16], ...
            'TooltipString', 'Start/stop animating between stored locations', ...
            'Callback', @onPlayPressed );
        
        % Remove some things we don't want from the toolbar and add a
        % toggle to the toolbar to hide the controls
        tb = findall( out.Window, 'Type', 'uitoolbar' );
        delete( findall( tb, 'Tag', 'Standard.FileOpen' ) );
        delete( findall( tb, 'Tag', 'Standard.NewFigure' ) );
        delete( findall( tb, 'Tag', 'Standard.EditPlot' ) );
        delete( findall( tb, 'Tag', 'Exploration.Brushing' ) );
        delete( findall( tb, 'Tag', 'Exploration.DataCursor' ) );
        delete( findall( tb, 'Tag', 'Exploration.Rotate' ) );
        delete( findall( tb, 'Tag', 'DataManager.Linking' ) );
        delete( findall( tb, 'Tag', 'Plottools.PlottoolsOn' ) );
        delete( findall( tb, 'Tag', 'Plottools.PlottoolsOff' ) );
        out.AnimateToggle = uitoggletool( ...
            'Parent', tb, ...
            'CData', readIcon( 'icon_play.png' ), ...
            'TooltipString', 'Start/stop animating between stored locations', ...
            'State', 'off', ...
            'Separator', 'on', ...
            'ClickedCallback', @onPlayToolbarPressed );
        out.ShowControlsToggle = uitoggletool( ...
            'Parent', tb, ...
            'CData', readIcon( 'icon_mandelControls.png' ), ...
            'TooltipString', 'Show/hide the Mandelbrot control panel', ...
            'State', 'on', ...
            'ClickedCallback', @onControlsTogglePressed );
        
        
        % Add listeners so that we can redraw when the axes are moved
        axHandle = handle( out.MainAxes );
        out.Listeners = [
            handle.listener( axHandle, findprop( axHandle, 'YLim' ), 'PropertyPostSet', @onLimitsChanged )
            ]; %#ok<NBRAK>
        % Also redraw if resized
        set( out.Window, 'ResizeFcn', @onFigureResize, ...
            'CloseRequestFcn', @onFigureClose );
    end % createGUI

    function onLimitsChanged( ~, ~ )
        redraw();
    end % onLimitsChanged

    function onFigureResize( ~, ~ )
        % Change the axes limits to exactly fit the figure
        pos = get( gui.Window, 'Position' );
        xlim = get( gui.MainAxes, 'XLim' );
        ylim = get( gui.MainAxes, 'YLim' );
        delta_ylim = ( diff( xlim )*pos(4)/pos(3) - diff( ylim ) ) / 2;
        data.WindowPixelSize = pos(3:4);
        % Set the YLim to give the correct aspect. This will trigger a
        % redraw
        set( gui.MainAxes, 'YLim', ylim + delta_ylim*[-1 1] );
    end % onFigureResize

    function onFigureClose( ~, ~ )
        % Clear up
        data.IsAnimating = false;
        if data.WriteVideo
            close( data.VideoWriter );
        end
        delete( gui.Window );
    end % onFigureClose

    function onCalculationMethodChanged( ~, ~ )
        idx = get( gui.MethodSelector, 'Value' );
        data.SelectedCalculationMethod = data.CalculationMethods{idx};
        redraw();
    end % onCalculationMethodChanged

    function onAddPressed( ~, ~ )
        disp('Add')
        fprintf( 'Adding location [%1.15f, %1.15f], %1.15f\n', ...
            data.XLim, data.Y );
        idx = numel( data.LocationList ) + 1;
        data.LocationList(idx).XLim = data.XLim;
        data.LocationList(idx).Y = data.Y;
        writeLocationList( data.LocationList );
        % Release the button
        set( gui.AddButton, 'Value', 0 );
    end % onAddPressed

    function onResetPressed( ~, ~ )
        disp('Reset')
        fprintf( 'Leaving location [%1.15f, %1.15f], %1.15f\n', ...
            data.XLim, data.Y );
        pos = get( gui.Window, 'Position' );
        aspect = pos(4)/pos(3);
        ylim = diff( data.OrigXLim ) * aspect / 2 * [-1 1];
        set( gui.MainAxes, 'XLim', data.OrigXLim, 'YLim', data.OrigY + ylim );
        set( gui.ResetButton, 'Value', 0 );
    end % onResetPressed

    function onPlayPressed( ~, ~ )
        disp('Play')
        if get( gui.AnimateButton, 'Value' )==1
            updateAnimationControls(true);
            while ishandle(gui.AnimateButton) && (get( gui.AnimateButton, 'Value' )==1)
                newXLim = data.LocationList(data.NextLocation).XLim;
                newY    = data.LocationList(data.NextLocation).Y;
                animatedMove( newXLim, newY );
                if numel( data.LocationList )>1
                    % Choose a random location
                    thisLocation = data.NextLocation;
                    while data.NextLocation == thisLocation
                        data.NextLocation = randi( numel( data.LocationList ), 1 );
                    end
                    fprintf( 'Next location: %d\n', data.NextLocation )
                else
                    % Only one location, so stop
                    data.IsAnimating = false;
                    updateAnimationControls( false );
                end
            end
        else
            data.IsAnimating = false;
            updateAnimationControls( false );
        end
    end % onPlayPressed

    function onPlayToolbarPressed( ~, evt )
        ison = strcmpi(get( gui.AnimateToggle, 'State' ), 'on');
        updateAnimationControls( ison );
        onPlayPressed(gui.AnimateButton, evt);
    end % onPlayToolbarPressed

    function updateAnimationControls( isAnimating )
        if isAnimating
            set( gui.AnimateButton, 'Value', 1 );
            set( gui.AnimateToggle, 'State', 'on' );
        else
            set( gui.AnimateButton, 'Value', 0 );
            set( gui.AnimateToggle, 'State', 'off' );
        end
        drawnow();
    end % updateAnimationControls

    function onControlsTogglePressed( ~, ~ )
        % Toggle the control panel on and off
        disp('Toggle controls')
        pos = get( gui.ControlPanel, 'Position' );
        if strcmpi( get( gui.ShowControlsToggle, 'State' ), 'off' )
            % Turn it off (move offscreen)
            pos(1) = -pos(3)-10;
        else
            % Turn it on (move onscreen)
            pos(1) = 10;
        end
        set( gui.ControlPanel, 'Position', pos );
        
    end

    function animatedMove( targetXLim, targetY )
        % Form a zoom path between the two
        data.IsAnimating = true;
        if isequal( data.XLim, targetXLim ) && isequal( data.Y, targetY )
            data.IsAnimating = false;
            return;
        end
        
        % Perform a zoom and translate arc
        maxNumSteps = 1000;
        distTravelled = sqrt( (mean( data.XLim ) - mean( targetXLim )).^2 ...
            + (data.Y - targetY).^2 );
        adjustRatio = exp( -10*linspace(-3,3,maxNumSteps).^2 );
        adjustRatio = adjustRatio - min(adjustRatio);
        ratio = cumsum( adjustRatio ); ratio = ratio / ratio(end);
        minXPath = interp1( [0,1], [data.XLim(1),targetXLim(1)], ratio );
        maxXPath = interp1( [0,1], [data.XLim(2),targetXLim(2)], ratio );
        maxXRange = max( maxXPath - maxXPath );
        xlimAdjust = max(0, 0.3*distTravelled - maxXRange);
        minXPath = minXPath - xlimAdjust*adjustRatio;
        maxXPath = maxXPath + xlimAdjust*adjustRatio;
        
        % Cull the ends if there's negligable motion. This helps to keep
        % things smooth but without long periods of no apparant motion.
        tolerance = 0.001;
        xRange = maxXPath - minXPath;
        firstGood = find( (xRange > (1+tolerance)*xRange(1)) | (xRange < (1-tolerance)*xRange(1)), 1, 'first' );
        if ~isempty( firstGood ) && firstGood > 2
            toCull = 2:firstGood-1;
        else
            toCull = [];
        end
        lastGood = find( (xRange > (1+tolerance)*xRange(end)) | (xRange < (1-tolerance)*xRange(end)), 1, 'last' );
        if ~isempty( lastGood ) && lastGood < numel(xRange)-1
            toCull = [toCull, lastGood:numel(xRange)-1];
        end
        ratio(toCull) = [];
        minXPath(toCull) = [];
        maxXPath(toCull) = [];
        xRange(toCull) = [];
        
        if ~isempty( ratio )
            % Work out the aspect ratio
            pos = get( gui.Window, 'Position' );
            aspect = pos(4)/pos(3);
            
            YPath = interp1( [0,1], [data.Y,targetY], ratio );
            heightPath = aspect*xRange;
            minYPath = YPath - 0.5*heightPath;
            maxYPath = YPath + 0.5*heightPath;
            for ii=1:numel(ratio)
                % Setting the limits will cause a redraw
                set( gui.MainAxes, ...
                    'XLim', [minXPath(ii),maxXPath(ii)], ...
                    'YLim', [minYPath(ii),maxYPath(ii)] );
                if data.IsAnimating == false
                    break;
                end
            end
        end
        data.IsAnimating = false;
        % Do a final redraw at full res
        redraw();
    end % animatedMove

    function kernel = loadKernel()
        thisDir = fileparts( mfilename( 'fullpath' ) );
        baseName = 'mandelbrotViewerProcessElement';
        data.CUDAFile = fullfile( thisDir, [baseName,'.cu'] );
        ptxname = [baseName,'.',parallel.gpu.ptxext];
        data.PTXFile = fullfile( thisDir, ptxname );
        if exist( data.PTXFile, 'file' ) ~= 2
            close( gui.Window );
            error( 'mandelbrotViewer:MissingPTX', 'Could not find ''%s''. Please use NVCC to compile it', ptxname );
        end
        kernel = parallel.gpu.CUDAKernel( data.PTXFile, data.CUDAFile );
    end % loadKernel

    function redraw()
        % Protect against the window closing
        if ~ishandle(gui.MainAxes)
            return;
        end
        % To work out what to draw and at what resolution we need the axis
        % limits and pixel counts.
        xlim = get(gui.MainAxes,'XLim');
        ylim = get(gui.MainAxes,'YLim');
        data.XLim = xlim;
        data.Y = mean( ylim );
        imWidth = data.WindowPixelSize(1);
        imHeight = data.WindowPixelSize(2);
        if data.IsAnimating && (imWidth*imHeight>600000)
            % To speed up animations with large windows, subsample by 2
            imWidth = round(imWidth/2);
            imHeight = round(imHeight/2);
        end
        
        zoomLevel = imWidth / diff( xlim );
        maxIterations = min( data.MaxIterations, 200 + 0.1*sqrt(zoomLevel) );
        
        % Call the computation
        t = tic;
        switch( data.SelectedCalculationMethod )
            case 'CUDAKernel'
                logCount = computeMandelbrotCUDAKernel( xlim, imWidth, ...
                    ylim, imHeight, ...
                    maxIterations );
                
            case 'CPU'
                logCount = computeMandelbrotCPU( xlim, imWidth, ...
                    ylim, imHeight, ...
                    maxIterations );
                
            case 'GPU (simple)'
                logCount = computeMandelbrotGPU( xlim, imWidth, ...
                    ylim, imHeight, ...
                    maxIterations );
                
            case 'GPU ArrayFun'
                logCount = computeMandelbrotArrayFun( xlim, imWidth, ...
                    ylim, imHeight, ...
                    maxIterations );
                
            otherwise
                error( 'mandelbrotViewer:BadMethod', 'Unrecognised calculation method ''%s''', data.SelectedCalculationMethod );
        end
        
        computeTime = toc(t);
        
        minCount = min( logCount(:) );
        logCount = (logCount - minCount) ./ (log(maxIterations+1)-minCount);
        % Guard against a closed window
        if ~ishandle( gui.Image )
            return;
        end
        set( gui.Image, ...
            'XData', xlim, ...
            'YData', ylim, ...
            'CData', logCount );
        if data.ControlsVisible
            set( gui.ComputeText, 'String', sprintf( 'Computed in %dms', round(1000*computeTime) ) )
            
            % Capture the current time for frame-rate calculations
            thisFrameTime = now();
            framerate = 1 / (86400*(thisFrameTime - data.LastFrameTime)); % convert days to seconds
            set( gui.FrameRateText, 'String', sprintf( 'Displaying at %dfps', round(framerate) ) )
            data.LastFrameTime = thisFrameTime;
            
            % Force a redraw
            drawnow();
        end
        
        % Capture!
        if data.WriteVideo
            t0 = now();
            currFrame = getframe( gui.Window );
            writeVideo( data.VideoWriter, currFrame );
            % Also reset the frame time to exclude the video writing
            delta_t = now() - t0;
            data.LastFrameTime = data.LastFrameTime + delta_t;
        end
        
        
    end % redraw


    function locations = readLocationList()
        fid = fopen( 'locations.csv', 'rt' );
        if fid<0
            close( gui.Window );
            error( 'mandelbrotViewer:BadLocationRead', 'Could not open location list for reading: ''locations.csv''' );
        end
        
        locData = textscan( fid, '%f,%f,%f' );
        N = size( locData{1}, 1 );
        if N<1
            close( gui.Window );
            error( 'mandelbrotViewer:EmptyLocationFile', 'No locations found in: ''locations.csv''' );            
        end
        locations = struct( ...
            'XLim', cell( N, 1 ), ...
            'Y', cell( N, 1 ) );
        for ii=1:N
            locations(ii).XLim = [locData{1}(ii), locData{2}(ii)];
            locations(ii).Y = locData{3}(ii);
        end
        
        fclose( fid );
    end % readLocationList

    function writeLocationList( locations )
        fid = fopen( 'locations.csv', 'wt' );
        if fid<0
            error( 'mandelbrotViewer:BadLocationWrite', 'Could not open location list for writing: ''locations.csv''' );
        end
        N = numel( locations );
        for ii=1:N
            fprintf( fid, '%1.15f,%1.15f,%1.15f\n', ...
                locations(ii).XLim(1), ...
                locations(ii).XLim(2), ...
                locations(ii).Y );
        end
        
        fclose( fid );
    end % writeLocationList

    function cdata = readIcon( filename )
        [cdata,~,alpha] = imread( filename );
        idx = find( ~alpha );
        page = size(cdata,1)*size(cdata,2);
        cdata = double( cdata ) / 255;
        cdata(idx) = nan;
        cdata(idx+page) = nan;
        cdata(idx+2*page) = nan;
    end % readIcon

    function matlabVersionCheck()
        % R2011a is v7.12
        majorMinor = sscanf( version, '%d.%d' );
        if (majorMinor(1)<7) || (majorMinor(1)==7 && majorMinor(2)<13)
            error( 'mandelbrotViewer:MATLABTooOld', 'mandelbrotViewer requires MATLAB R2011b or above.' );
        end
    end % matlabVersionCheck

    function gpuCheck()
        try
            d = gpuDevice();
        catch err
            error( 'mandelbrotViewer:NoGPU', 'mandelbrotViewer requires a GPU and none appear to be availble. Type "gpuDevice" for more information.' );
        end
        if ~d.DeviceSupported
            error( 'mandelbrotViewer:GPUNotSupported', 'The selected GPU is not supported. Type "gpuDevice" for more information.' );
        end
    end % matlabVersionCheck

end % mandelbrotViewer