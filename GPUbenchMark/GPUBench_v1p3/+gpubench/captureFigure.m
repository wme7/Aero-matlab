function captureFigure( figh, filename, trim )
%captureFigure: use print to capture a figure with over-sampling
%
%   gpubench.captureFigure(figh,filename,trim)

%   Copyright 2011 The MathWorks, Inc.

if nargin<3
    trim = true;
end

% Have to hard-code the DPI since Windows doesn't report it correctly
% dpi = get(0,'ScreenPixelsPerInch');
%dpi = 192;
dpi = 160;

% Copy the position to the paper position to get the size right
pos = get( figh, 'Position' );
set( figh, 'PaperUnits', 'Points', 'PaperPosition', pos );

% Print the image to file at screen resolution
print( figh, '-zbuffer', sprintf('-r%d',round(0.5*dpi)), '-dpng', filename );

% Load and trim the top whitespace
gpubench.trimImage( filename, trim );
