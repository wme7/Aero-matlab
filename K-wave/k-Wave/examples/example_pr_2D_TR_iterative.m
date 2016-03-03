% 2D Iterative Image Improvement Using Time Reversal Example
%
% This example demonstrates how photoacoustic image reconstruction may be
% improved iteratively. First, the sensor data is simulated using
% kspaceFirstOrder2D. Then, in the first reconstruction, an image is formed
% from data recorded on a line array using time reversal. In the second, an
% image is formed from data recorded on an L-shaped sensor array (also
% using time reversal). In the final stage, this image is improved
% iteratively. This example builds on the 2D Time Reversal Reconstruction
% For A Line Sensor Example.
%
% author: Ben Cox
% date: 22nd January 2012
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

clear all;

% =========================================================================
% SET UP THE SIMULATION
% =========================================================================

% create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 128 - 2*PML_size;  % number of grid points in the x (row) direction
Ny = 256 - 2*PML_size;  % number of grid points in the y (column) direction
dx = 0.1e-3;            % grid point spacing in the x direction [m]
dy = 0.1e-3;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% load an image for the initial pressure distribution
p0_image = loadImage('EXAMPLE_k-Wave.png');

% make it binary
p0_image = (p0_image>0);

% smooth and scale the initial pressure distribution
p0_magnitude = 2;
p0 = p0_magnitude*smooth(kgrid, p0_image, true);

% assign to the source structure
source.p0 = p0;

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, 'PlotPML', false};

% =========================================================================
% DEFINE SENSOR INDICES FOR LATER USE
% =========================================================================

% define a four-sided box array and find the indices of the sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(1,:) = 1;
sensor.mask(:,1) = 1;
sensor.mask(end,:) = 1;
sensor.mask(:,end) = 1;
sensor_indices = find(sensor.mask==1);

% find the indices along two sides, with respect to sensor_indices
sensor.mask = zeros(Nx, Ny);
sensor.mask(end,:) = 1;
sensor.mask(:,end) = 1;
sensor_box_indices1 = find(sensor.mask(sensor_indices)==1);

% define an L-shaped array used for the measurements, find the sensor
% indices, and their indices with respect to sensor_indices 
sensor.mask = zeros(Nx, Ny);
sensor.mask(1,:) = 1;
sensor.mask(:,1) = 1;
original_sensor_indices = find(sensor.mask == 1);
sensor_box_indices2 = find(sensor.mask(sensor_indices)==1);

% =========================================================================
% RUN THE SIMULATION AND RECORD THE DATA OVER THE L-SHAPED ARRAY
% =========================================================================

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% RECONSTRUCT AN IMAGE USING TIME REVERSAL AND LINE ARRAY DATA
% =========================================================================

% define a line array of sensors
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data((sensor.mask(original_sensor_indices) == 1),:);

% set the initial pressure to be zero
source.p0 = 0;

% run the time reversal reconstruction
p0_line = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% RECONSTRUCT AN IMAGE USING TIME REVERSAL AND DATA FROM THE L-SHAPED ARRAY
% =========================================================================

% extend the array to the original L-shaped form
sensor.mask(:, 1) = 1;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_L = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% ITERATIVELY IMPROVE THE IMAGE
% =========================================================================

% begin with the image already recovered 
p0_iterated = p0_L;

Niterations = 5;
for loop = 1:Niterations

    % apply a non-negativity condition 
    source.p0 = p0_iterated.*(p0_iterated>=0);
    
    % set the sensor mask to be just the missing sides of the box
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(end, :) = 1;
    sensor.mask(:,end) = 1;

    % remove the time-reversed field so the forward model will run
    sensor = rmfield(sensor,'time_reversal_boundary_data');

    % estimate the data on the missing sides of the box
    estimated_sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    
    % define the full box array
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(1,:) = 1;
    sensor.mask(:,1) = 1;
    sensor.mask(end,:) = 1;
    sensor.mask(:,end) = 1;
    
    % combine the measured data with the newly estimated data
    combined_sensor_data = zeros(sum(sensor.mask(:)==1), length(kgrid.t_array));
    combined_sensor_data(sensor_box_indices1,:) = estimated_sensor_data;
    combined_sensor_data(sensor_box_indices2,:) = sensor_data;
    
    % assign the time reversal data
    sensor.time_reversal_boundary_data = combined_sensor_data;    
    
    % reset the initial pressure
    source.p0 = 0;

    % run the time reversal reconstruction
    p0_iterated = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure
subplot(2,1,1)
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p0, [-p0_magnitude p0_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('true image')
axis image;
colorbar;

% plot the initial pressure reconstructed using a line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1,:) = 1;
subplot(2,1,2)
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p0_line + sensor.mask*p0_magnitude, [-p0_magnitude p0_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('line sensor')
axis image;
colorbar;

% plot the initial pressure reconstructed using an L-shaped sensor
sensor.mask(:,1) = 1;
figure
subplot(2,1,1)
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p0_L + sensor.mask*p0_magnitude, [-p0_magnitude p0_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('L-shaped sensor')
axis image;
colorbar;

% plot the initial pressure reconstructed iteratively
subplot(2,1,2)
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, p0_iterated + sensor.mask*p0_magnitude, [-p0_magnitude p0_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
title('Iterated')
axis image;
colorbar;

% plot a profile for comparison
figure;
plotline = 45;
plot(kgrid.y_vec*1e3, p0(plotline, :), 'k-', ...
    kgrid.y_vec*1e3, p0_line(plotline, :), 'r-', ...
    kgrid.y_vec*1e3, p0_L(plotline, :), 'm-', ...
    kgrid.y_vec*1e3, p0_iterated(plotline, :), 'b--');
xlabel('y-position [mm]');
ylabel('Pressure');
legend('Initial Pressure', 'Line', 'L-shaped', 'Iterated');
axis tight;
set(gca, 'YLim', [0 2.6]);
