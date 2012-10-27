function [logCount,t] = mandelbrotViewerProcessElement( x0, y0, escapeRadius2, maxIterations )
% Evaluate the Mandelbrot function for a single element

%   Copyright 2010-2011 The Mathworks, Inc.

t = 1;
z0 = complex( x0, y0 );
z = z0;
count = 0;
while count <= maxIterations && (z*conj(z) <= escapeRadius2)
    z = z*z + z0;
    count = count + 1;
end
magZ2 = max(real(z).^2 + imag(z).^2,escapeRadius2);
logCount = log( count + 1 - log( log( magZ2 ) / 2 ) / log(2) );
