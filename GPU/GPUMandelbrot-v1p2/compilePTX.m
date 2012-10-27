function compilePTX()
%compilePTX  try to recompile the PTX for the Mandelbrot kernel
%
%   compilePTX() attempts to call NVCC to recompile the Mandelbrot "process
%   element" kernel for the current platform. The PTX is produced for 1.3
%   cards to allow use on any supported device. If you have a 2.0 (Fermi)
%   or 3.0 (Kepler) device, you may see some performance benefit to
%   chnageing the architecture flag below to match.

%   Copyright 2012 The Mathworks, Inc.

ptxarch = 'sm_13';
kernel = 'mandelbrotViewerProcessElement';

flags = sprintf( '-arch=%s', ptxarch );
if ismac
    % On Mac we must force the use of 64-bit pointers otherwise the host
    % uses 64-bit and the device 32!
    flags = [flags, ' -m 64'];
end

cmd = sprintf( 'nvcc -ptx %s.cu %s -o %s.%s', ...
    kernel, flags, kernel, parallel.gpu.ptxext );
fprintf('Running command: %s\n', cmd );
result = system( cmd );
if (result ~= 0)
    fprintf( 'Failed with error %d.\n', result )
end
end
    