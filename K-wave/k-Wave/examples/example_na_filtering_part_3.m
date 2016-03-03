% Filtering A Delta Function Input Signal Example Part 3
%
% This example illustrates how to temporally filter an input signal.
%
% author: Bradley Treeby
% date: 19th January 2010
% last update: 24th August 2014
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

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

clear all

% modify this parameter to run the different examples
example_number = 1;
% 1: default causal filter
% 2: zero phase filter
% 3: causal filter with modified input settings

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 256;           % number of grid points in the x (row) direction
dx = 10e-3/Nx;      % grid point spacing in the x direction [m]
kgrid = makeGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create a time array
num_time_steps = 1024;
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
kgrid.t_array = 0:dt:dt*(num_time_steps - 1);

% define a delta function input pulse
temporal_offset = 100;      % [time steps]
source_magnitude = 2;       % [au]
source_func = zeros(size(kgrid.t_array));
source_func(temporal_offset) = source_magnitude;

switch example_number
    case 1
        % filter the input signal
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func);
    case 2
        % filter the input signal
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func, 'ZeroPhase', true);
    case 3
        % filter the input signal
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func, 'PPW', 4, 'TransitionWidth', 0.05);
end

% compute the amplitude spectra of the original and filtered input signals
[f, source_func_as] = spect(source_func, 1/dt);
[f, source_func_filtered_as] = spect(source_func_filtered, 1/dt);

% extract the maximum frequency supported by the grid (two points per
% wavelength)
f_max = kgrid.k_max * min(medium.sound_speed(:)) / (2*pi);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot original and filtered input signals
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
subplot(2, 1, 1), plot(kgrid.t_array*scale, source_func, 'k-', kgrid.t_array*scale, source_func_filtered, 'b-');
xlabel(['Time [' prefix 's]']);
ylabel('Pressure [au]');
legend('original input', 'filtered input');
axis tight;
set(gca, 'XLim', [0 3]);

% plot the amplitude spectra
[f_sc, scale, prefix] = scaleSI(max(f));
subplot(2, 1, 2), plot(f*scale, source_func_as, 'k-', f*scale, source_func_filtered_as, 'b-');
xlabel(['Frequency [' prefix 'Hz]']);
ylabel('Amplitude [au]');

% plot the maximum frequency supported by the grid
ylim = get(gca, 'YLim');
hold on;
subplot(2, 1, 2), line([f_max*scale, f_max*scale], [0, ylim(2)], 'LineStyle','--', 'Color', 'k');
set(gca, 'XLim', [0 50]);