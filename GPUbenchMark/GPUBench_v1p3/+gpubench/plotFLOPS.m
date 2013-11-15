function plotFLOPS( data, color, linewidth )
%plotFLOPS  plot some gpuBench data as floating-point operations per second
%
%   gpubench.plotFLOPS(data,color,linewidth)

%   Copyright 2011 The MathWorks, Inc.

semilogx( data.Sizes, 1e-9 * data.NumOps ./ data.Times, ...
    'Color', color, ...
    'Marker', '.', ...
    'MarkerSize', 10+6*linewidth, ...
    'linewidth', linewidth )
xlabel( 'Number of elements' )
ylabel( sprintf('GFLOPS\n(higher is better)') )
grid( 'on' )
hold( 'on' )
