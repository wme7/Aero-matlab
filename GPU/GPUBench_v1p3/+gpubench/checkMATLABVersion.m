function checkMATLABVersion()
%checkMATLABVersion  error if the MATLAB version is not suitable

%   Copyright 2011 The MathWorks, Inc.

v = ver('MATLAB');
v = regexp( v.Version, '(?<major>\d+)\.(?<minor>\d+)', 'names' );
v.major = str2double( v.major );
v.minor = str2double( v.minor );
if (v.major < 7) || ((v.major==7) && (v.minor<13))
    error( 'GPUBench:MATLABVersion', 'GPUBench requires MATLAB version 7.13 (R2011b) or higher' );
end
