function data = timeShift(data, dt, forward)
%TIMESHIFT   Shift time series to and from staggered temporal grid.
%
% DESCRIPTION:
%       timeShift shifts an input time series (or 2D array of time series)
%       to and from a staggered temporal grid using Fourier interpolation.
%       If the input data is defined as a 2D array of time series, the
%       temporal shift is performed along the second dimension assuming
%       that the array is indexed as (n, t).
%
%       This function can be used to shift the acoustic particle velocity
%       recorded by the first-order simulation functions to the regular
%       (non-staggered) temporal grid.
%
%       Example:
%
%           dt      = pi/5;
%           t       = 0:dt:(2*pi - dt);
%           y       = cos(t);
%           y_shift = timeShift(y, dt);
%           y_ref   = cos(t + dt/2);
%           plot(t, y, 'k-s', t, y_shift, 'b-s', t, y_ref, 'rx');
%
% USAGE:
%       data = timeShift(data, dt)
%       data = timeShift(data, dt, true)
%       data = timeShift(data, dt, false)
%
% INPUTS:
%       data        - 1D or 2D input data, for 2D input data, the time
%                     shift is performed along the second dimension
%       dt          - time step [s]
%
% OPTIONAL INPUTS:
%       forward     - Boolean controlling whether the time series is
%                     shifted forward or backward in time [default = true]
%
% OUTPUTS:
%       data        - shifted time series
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 5th September 2013
%       last update - 5th September 2013
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox
%
% See also gradientSpect, kspaceFirstOrder1D, kspaceFirstOrder2D,
% kspaceFirstOrder3D, pstdElastic2D, pstdElastic3D

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% check the size of the input
if numDim(data) > 2
    error('input data must be 1 or 2 dimensional');
end

% check for optional "forward" input
if nargin < 3
    forward = true;
end

% define the discretisation of the temporal dimension such that there is
% always a DC component 
Nt = size(data, 2);
if rem(Nt, 2) == 0
    % grid dimension has an even number of points
    w_vec = ((-Nt/2:Nt/2-1)/Nt);
else
    % grid dimension has an odd number of points
    w_vec = ((-(Nt-1)/2:(Nt-1)/2)/Nt);
end

% force middle value to be zero in case 1/Nt is a recurring number and the
% series doesn't give exactly zero 
w_vec(floor(Nt/2) + 1) = 0;

% define the temporal frequency vector
w_vec = (2*pi/dt) .* w_vec;  

% create the temporal shift operator
if forward
    t_shift_pos = ifftshift( exp(1i*w_vec*dt/2) );
else
    t_shift_pos = ifftshift( exp(-1i*w_vec*dt/2) );
end

% shift the input to the regular (non-staggered) temporal grid
switch numDim(data)
    case 1
        % shift along non-singleton dimension
        data = real(ifft(bsxfun(@times, t_shift_pos, fft(data))));
    case 2
        % shift along second dimension assuming data is indexed as (n, t) 
        data = real(ifft(bsxfun(@times, t_shift_pos, fft(data, [], 2)), [], 2));
end
