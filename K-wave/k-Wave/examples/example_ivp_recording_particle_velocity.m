% Recording Particle Velocity Example
%
% This example demonstrates how to record the particle velocity using a
% Cartesian or binary sensor mask. It builds on the Homogeneous Propagation
% Medium and Heterogeneous Propagation Medium examples.  
%
% author: Bradley Treeby
% date: 1st November 2010
% last update: 20th September 2012
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
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.density = 1000;      % [kg/m^3]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5; 

% create time array
t_end = 6e-6;
kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [au]
disc_x_pos = Nx/2;  % [grid points]
disc_y_pos = Ny/2;  % [grid points]
disc_radius = 5;    % [grid points]
source.p0 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define four sensor points centered about source.p0
sensor_radius = 40; % [grid points]
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2 - sensor_radius, Ny/2) = 1;
sensor.mask(Nx/2, Ny/2 + sensor_radius) = 1;
sensor.mask(Nx/2, Ny/2 - sensor_radius) = 1;

% set the acoustic variables that are recorded
sensor.record = {'p', 'u'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, source.p0 + sensor.mask, [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the simulated sensor data
figure;
[t, t_sc, t_prefix] = scaleSI(kgrid.t_array(end));
mx = 5e-7;
for sensor_num = 1:4
    % plot the pressure
    subplot(4, 3, 3*sensor_num - 2), plot(t_sc*kgrid.t_array, sensor_data.p(sensor_num, :), 'k-');
    set(gca, 'YLim', [-0.75, 0.75], 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('p');    

    % plot the particle velocity ux
    subplot(4, 3, 3*sensor_num - 1), plot(t_sc*kgrid.t_array, sensor_data.ux(sensor_num, :), 'k-');
    set(gca, 'YLim', [-mx, mx], 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('ux'); 

    % plot the particle velocity uz
    subplot(4, 3, 3*sensor_num), plot(t_sc*kgrid.t_array, sensor_data.uy(sensor_num, :), 'k-');
    set(gca, 'YLim', [-mx, mx], 'XLim', [0 t_end*t_sc]);
    xlabel(['time [' t_prefix 's]']);
    ylabel('uy'); 
end
