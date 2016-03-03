% Shear Waves And Critical Angle Reflection Example
%
% This example illustrates Snell's law for elastic media using a weakly
% focused ultrasound transducer incident on a soft-tissue / bone interface.
% It builds on the Explosive Source In A Layered Medium and Snell's Law And
% Critical Angle Reflection examples.
%
% author: Bradley Treeby
% date: 14th Nov 2013
% last update: 15th February 2014
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
% SIMULATION PARAMETERS
% =========================================================================

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

% create the computational grid
PML_size = 10;               % size of the PML in grid points
Nx = 128*scale - 2*PML_size; % number of grid points in the x direction
Ny = 192*scale - 2*PML_size; % number of grid points in the y direction
dx = 0.1e-3/scale;           % grid point spacing in the x direction [m]
dy = 0.1e-3/scale;           % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the medium properties
cp1 = 1540;         % compressional wave speed [m/s]
cs1 = 0;            % shear wave speed [m/s]
rho1 = 1000;        % density [kg/m^3]
alpha0_p1 = 0.1;    % compressional wave absorption [dB/(MHz^2 cm)]
alpha0_s1 = 0.1;    % shear wave absorption [dB/(MHz^2 cm)]

cp2 = 3000;         % compressional wave speed [m/s]
cs2 = 1400;         % shear wave speed [m/s]
rho2 = 1850;        % density [kg/m^3]
alpha0_p2 = 1;      % compressional wave absorption [dB/(MHz^2 cm)]
alpha0_s2 = 1;      % shear wave absorption [dB/(MHz^2 cm)]

% create the time array
cfl   = 0.1;
t_end = 12e-6;
kgrid.t_array= makeTime(kgrid, cp1, cfl, t_end);

% define position of heterogeneous slab
slab = zeros(Nx, Ny);
slab(Nx/2:Nx, :) = 1;

% define the source properties
source_freq = 2e6;              % [Hz]
source_strength = 1e6;          % [Pa]
source_cycles = 3;              % number of tone burst cycles
source_focus_dist = 5*scale;    % position of focus inside slab [grid points]
source_slab_dist = 5*scale;     % distance between source and slab [grid points]
source_mask = makeCircle(Nx, Ny, Nx/2 + source_focus_dist, Ny/2, 50*scale, pi/3);
source_mask(Nx/2 - source_slab_dist:end, :) = 0;

% define the sensor to record the maximum particle velocity everywhere
sensor.record = {'u_max_all'};

% set the input arguments
input_args = {'PMLSize', PML_size, 'PMLAlpha', 2, 'PlotPML', false, ...
    'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
    'DisplayMask', 'off', 'DataCast', 'single'};

% =========================================================================
% FLUID SIMULATION
% =========================================================================

% assign the medium properties
medium.sound_speed            = cp1*ones(Nx, Ny);
medium.sound_speed(slab == 1) = cp2;
medium.density                = rho1*ones(Nx, Ny);
medium.density(slab == 1)     = rho2;
medium.alpha_coeff            = alpha0_p1*ones(Nx, Ny);
medium.alpha_coeff(slab == 1) = alpha0_p2;
medium.alpha_power            = 2;

% assign the source
source.p_mask = source_mask;
source.p = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);

% run the fluid simulation
sensor_data_fluid = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% ELASTIC SIMULATION
% =========================================================================

% define the medium properties
clear medium
medium.sound_speed_compression            = cp1*ones(Nx, Ny);
medium.sound_speed_compression(slab == 1) = cp2;
medium.sound_speed_shear                  = cs1*ones(Nx, Ny);
medium.sound_speed_shear(slab == 1)       = cs2;
medium.density                            = rho1*ones(Nx, Ny);
medium.density(slab == 1)                 = rho2;
medium.alpha_coeff_compression            = alpha0_p1*ones(Nx, Ny);
medium.alpha_coeff_compression(slab == 1) = alpha0_p2;
medium.alpha_coeff_shear                  = alpha0_s1*ones(Nx, Ny);
medium.alpha_coeff_shear(slab == 1)       = alpha0_s2;

% assign the source
clear source
source.s_mask = source_mask;
source.sxx = -source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);
source.syy = source.sxx;

% run the elastic simulation
sensor_data_elastic = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% define plot vector
x_vec = kgrid.x_vec(1 + PML_size:end - PML_size)*1e3;
y_vec = kgrid.y_vec(1 + PML_size:end - PML_size)*1e3;

% calculate square of velocity magnitude
u_e = sensor_data_elastic.ux_max_all.^2 + sensor_data_elastic.uy_max_all.^2;
u_f = sensor_data_fluid.ux_max_all.^2 + sensor_data_fluid.uy_max_all.^2;

% plot layout
figure;
imagesc(y_vec, x_vec, double(source_mask | slab));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colormap(flipud(gray));

% plot beam patterns
figure;
subplot(2, 1, 1);
imagesc(y_vec, x_vec, 20*log10(u_f./max(u_f(:))));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colorbar;
caxis([-50, 0]);
title('Fluid Model');

subplot(2, 1, 2);
imagesc(y_vec, x_vec, 20*log10(u_e./max(u_e(:))));
xlabel('y [mm]');
ylabel('x [mm]');
axis image;
colorbar;
caxis([-50, 0]);
title('Elastic Model');
colormap(jet(256));