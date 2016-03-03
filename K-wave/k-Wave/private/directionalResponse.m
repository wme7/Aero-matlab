function p_sensor = directionalResponse(kgrid, sensor, sensor_mask_index, p_k)
%DIRECTIONALRESPONSE    Apply directivity to the sensor measurements.
% 
% DESCRIPTION:
%       directionalResponse takes a pressure field in k-space, p_k, 
%       multiplies by a directivity function, and outputs the grid points
%       indicated in sensor.mask.
%
% USAGE:
%       p_sensor = directionalResponse(kgrid, sensor, p_k)
%
% INPUTS:
%       kgrid       - k-Wave grid structure returned by makeGrid
%
%       sensor      - k-Wave sensor structure containing the following
%                     fields:
%
%           .directivity_angle   
%                   - a matrix with the same structure as
%                     sensor.mask that allocates a directivity angle to
%                     each sensor element as defined in sensor.mask. The
%                     angles are in radians:
%                     0 = max sensitivity in y direction (up/down)
%                     pi/2 or -pi/2 = max sensitivity in x direction (left/right)
%
%       	.directivity_unique_angles
%                   - list of the unique directivity angles returned by
%                     sensor.directivity_unique_angles =
%                     unique(sensor.directivity_angle(sensor.mask == 1)); 
%
%       	.directivity_pattern 
%                   - a text string with currently only one
%                     option. 'pressure' indicates that the directional
%                     response should be of the kind due to spatial
%                     averaging over a sensor surface, so a sinc function
%                     in 2D.
%
%       	.directivity_size    
%                   - the directivity pattern used is what
%                     would be the directivity if the sensor were this
%                     length (width). The larger this is the more
%                     directional the response.
%
%       	.directivity_wavenumbers
%                   - this is set to [kgrid.ky(:)'; kgrid.kx(:)']. It is
%                     precomputed to allow data casting, as kgrid.kx (etc)
%                     are computed on the fly.
%
%       sensor_mask_index
%                   - indices of the active sensor elements in sensor.mask
%
%       p_k         - the acoustic pressure field in the k-space domain.
% 
% Currently works for binary sensor_mask, but not when sensor_mask is given
% as Cartesian coordinates. Also, currently works only in 2D.
% 
% ABOUT:
%       author: Ben Cox and Bradley Treeby
%       date: 21st January 2010
%       last update: 25th August 2014
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

% sensor mask indices
Ns = length(sensor_mask_index);

% pre-allocate p_sensor vector
p_sensor = zeros(Ns, 1, class(sensor.directivity_angle));

% loop over number of unique angles
for loop = 1:length(sensor.directivity_unique_angles);
   
   % get current angle
   theta = sensor.directivity_unique_angles(loop);

   % find which of the sensors have this directivity
   indices = find(sensor.directivity_angle(sensor_mask_index) == theta);

   switch sensor.directivity_pattern
      case 'pressure'
         % calculate magnitude of component of wavenumber along sensor face
         k_tangent = reshape([cos(theta) -sin(theta)]*sensor.directivity_wavenumbers, kgrid.Nx, kgrid.Ny);
         directionality = fftshift(sinc(k_tangent*sensor.directivity_size/2));         
      case 'gradient'
         % calculate magnitude of component of wavenumber normal to the sensor face          
         k_normal = reshape([sin(theta), cos(theta)]*[kgrid.ky(:)'; kgrid.kx(:)'], kgrid.Nx, kgrid.Ny);
         temp = k_normal./kgrid.k;
         temp(kgrid.k==0) = 0;
         directionality = fftshift(temp);         
      otherwise
         error('Unsupported directivity pattern')
   end
   
   % apply the directivity response to the pressure field (in k-space)
   p_directivity = real(ifft2(p_k.*directionality));
   
   % pick out the response at the sensor points
   p_sensor(indices) = p_directivity(sensor_mask_index(indices));
   
end