% Modelling Sensor Directivity in 3D Example 
%
% This example demonstrates how the sensitivity of a large single element
% detector varies with the angular position of a point-like source. It is a
% 3D version of the Modelling Sensor Directivity in 2D example.
%
% author: Ben Cox
% date: 29th October 2010
% last update: 23rd February 2012
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
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 100e-3/Nx;     % grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
dz = dx;            % grid point spacing in the z direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% create the time array
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
Nt = length(kgrid.t_array);

% define a large area detector 
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(Nx/2+1, (Ny/2-9):(Ny/2+11),(Nz/2-9):(Nz/2+11)) = 1;

% define equally spaced point sources lying on a circle centred at the
% centre of the detector face 
Nangles = 11;
circle = makeCartCircle(25*dx, Nangles, [0,0], pi);
circle = [circle; zeros(1,Nangles)];

% find the binary sensor mask most closely corresponding to the cartesian
% points coordinates from makeCartCircle
circle3D = cart2grid(kgrid,circle);

% find the indices of the sources in the binary source mask
source_positions = find(circle3D == 1);

% define a time varying sinusoidal source
source_freq = 0.25e6;
source_mag = 1;
source.p = source_mag*sin(2*pi*source_freq*kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% pre-allocate array for storing the output time series
single_element_data = zeros(Nt,length(source_positions));

% run a simulation for each of these sources to see the effect that the
% angle from the detector has on the measured signal
for source_loop = 1:length(source_positions)
    
    % select a point source
    source.p_mask = zeros(Nx,Ny,Nz);
    source.p_mask(source_positions(source_loop)) = 1;

    % create a display mask to display the transducer
    display_mask = source.p_mask + sensor.mask;

    % run the simulation
    input_args = {'PMLSize', 10, 'DisplayMask', display_mask, 'PlotScale', [-0.2 0.2], 'PlotFreq', 50,  'DataCast', 'single'};
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

    % average the data recorded for each grid point to simulate the
    % measured signal from a large aperture, single element, detector
    single_element_data(:,source_loop) = sum(sum(sensor_data,1),1);

end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot source points and sensor mask
voxelPlot(circle3D + sensor.mask);
view([-34 34])

% plot the time series recorded for each of the sources
figure;
plot(kgrid.t_array, single_element_data);
colormap(getColorMap);
xlabel('time [s]');
ylabel('pressure');
title('time series from each direction');

% calculate angle between source and centre of detector face
angles = atan((kgrid.y(source_positions))./kgrid.x(source_positions));

% plot the maximum amplitudes for each of the sources, showing that the
% detector sensitivity falls off at low angles as expected.
figure;
plot(angles,max(single_element_data),'o')
colormap(getColorMap);
xlabel('angle between source and centre of detector face');
ylabel('maximum detected pressure from each direction')
