% Sensor Element Directivity in 2D Example
%
% This example shows how to attribute a directional response to a
% single-element sensor, or individual elements of a multi-element sensor
% array. Directionality can be included without a separate function through
% explicit averaging, as shown in the examples Modelling Sensor Directivity
% in 2D and Focussed Detector in 2D, but the functionality described here
% allows greater flexibility. Note that directivity defined in this way is
% currently only supported in 2D. This example builds on the Homogeneous
% Propagation Medium and Using A Binary Sensor Mask examples.
%
% author: Ben Cox
% date: 20th January 2010
% last update: 22nd February 2012
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

% =========================================================================
% DEFINE THE GRID AND MEDIUM PROPERTIES
% =========================================================================

% create the computational grid
Nx = 64;        % number of grid points in the x (row) direction
Ny = 64;        % number of grid points in the y (column) direction
dx = 1e-3/Nx;   % grid point spacing in the x direction [m]
dy = dx;     	% grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create the time array, then halve it to make the example end sooner
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed);
kgrid.t_array = (0:0.7*length(kgrid.t_array))*dt;

% =========================================================================
% DEFINE THE DIRECTIONAL SENSOR ARRAY
% =========================================================================

% define a line of sensor points
sensor.mask = zeros(Nx,Ny);
sensor.mask(24,2:2:63) = 1;

% define the angle of max directivity for each sensor point:
%    0             = max sensitivity in x direction (up/down)
%    pi/2 or -pi/2 = max sensitivity in y direction (left/right)
sensor.directivity_angle = zeros(Nx,Ny);
sensor.directivity_angle(24,2:2:63) = (-1:1/15:1)*pi/2;

% define the directivity pattern
sensor.directivity_pattern = 'pressure';

% define the directivity size
sensor.directivity_size = 16*kgrid.dx;

% =========================================================================
% SIMULATION AND VISUALISATION FOR AN INITIAL VALUE PROBLEM
% =========================================================================

% define the initial pressure distribution
source.p0 = zeros(Nx,Ny);
source.p0(39:41, :) = 2;
 
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLAlpha', [2 0]);

% plot the largest value of the output for each sensor
figure;
plot(((-1:1/15:1)*pi/2), max(sensor_data, [], 2), 'o');
xlabel('sensor directivity angle (radians)')
ylabel('maxima of single-element sensors'' outputs')

% =========================================================================
% SIMULATION AND VISUALISATION FOR TIME-VARYING SOURCE PROBLEM
% =========================================================================

% define a time varying sinusoidal source (instead of an initial pressure)
source.p_mask = source.p0;
source.p0 = [];
source_freq = 12e6;
source_mag = 0.25;
source.p = source_mag*sin(2*pi*source_freq*kgrid.t_array);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLAlpha', [2 0]);

% plot the largest value of the output for each sensor
max_data = max(sensor_data(:, 150:end), [], 2);
figure, plot(((-1:1/15:1)*pi/2), max_data./max(max_data),'o')
xlabel('sensor directivity angle (radians)')
ylabel('maxima of single-element sensors'' outputs')
axis([-pi/2 pi/2 0 1.05])
figure, polar(((-1:1/15:1)'*pi/2), max_data./max(max_data))
