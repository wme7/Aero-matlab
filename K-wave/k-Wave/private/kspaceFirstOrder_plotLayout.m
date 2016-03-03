% DESCRIPTION:
%       subscript to plot the simulation layout
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 7th February 2012
%       last update - 22nd February 2012
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

% set the colormap to use in 2D and 3D
plot_layout_cmap = flipud(bone(256));

switch kgrid.dim
    case 1
        
        % set the spacing as a fraction of the plot range
        plot_spc = 0.05;
                
        % create figure window
        figure;
        
        % plot the initial pressure distribution or the source mask
        subplot(2, 2, 1);
        if isfield(source, 'p0')
            plot(kgrid.x_vec(x1:x2)*scale, double(source.p0(x1:x2)));
            axis tight;
            plot_ylim = get(gca, 'YLim');
            plot_rng = plot_ylim(2) - plot_ylim(1);
            set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);
            title('Initial Pressure');
        elseif p_source
            plot(kgrid.x_vec(x1:x2)*scale, double(source.p_mask(x1:x2)));
            axis tight;
            plot_ylim = get(gca, 'YLim');
            plot_rng = plot_ylim(2) - plot_ylim(1);
            set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);            
            title('Source Mask');
        elseif u_source
            plot(kgrid.x_vec(x1:x2)*scale, double(source.u_mask(x1:x2)));
            axis tight;
            plot_ylim = get(gca, 'YLim');
            plot_rng = plot_ylim(2) - plot_ylim(1);
            set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);            
            title('Source Mask');            
        end
        
        % plot the sensor mask
        if record.use_sensor
            subplot(2, 2, 2), bar(kgrid.x_vec(x1:x2)*scale, double(sensor.mask(x1:x2)), 'b');
            axis tight;
            plot_ylim = get(gca, 'YLim');
            plot_rng = plot_ylim(2) - plot_ylim(1);
            set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);
            title('Sensor Mask');
        end
        
        % plot the sound speed distribution
        subplot(2, 2, 3);
        if numel(c) == 1
            plot(kgrid.x_vec(x1:x2)*scale, double(c));
        else
            plot(kgrid.x_vec(x1:x2)*scale, double(c(x1:x2)));
        end
        axis tight;
        plot_ylim = get(gca, 'YLim');
        plot_rng = plot_ylim(2) - plot_ylim(1);
        set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);
        title('Sound Speed');
        
        % plot the density distribution
        subplot(2, 2, 4);
        if numel(rho0) == 1
            plot(kgrid.x_vec(x1:x2)*scale, double(rho0));
        else
            plot(kgrid.x_vec(x1:x2)*scale, double(rho0(x1:x2)));
        end
        axis tight;
        plot_ylim = get(gca, 'YLim');
        plot_rng = plot_ylim(2) - plot_ylim(1);
        set(gca, 'YLim', [plot_ylim(1) - plot_spc*plot_rng, plot_ylim(2) + plot_spc*plot_rng]);
        title('Density');
        
        % add axis label
        xlabel(['(All horizontal axes in ' prefix 'm)']);
               
    case 2
        
        % set the spacing as a fraction of the plot range
        plot_spc = 0.6;        
        
        % create figure window
        figure;        

        % plot the initial pressure distribution or the source mask
        if isfield(source, 'p0')
            subplot(2, 2, 1), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(source.p0(x1:x2, y1:y2)));
            axis image;
            colormap(plot_layout_cmap);
            title('Initial Pressure');            
        elseif p_source
            subplot(2, 2, 1), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(source.p_mask(x1:x2, y1:y2)));
            axis image;
            colormap(plot_layout_cmap);
            title('Source Mask');
        elseif u_source
            subplot(2, 2, 1), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(source.u_mask(x1:x2, y1:y2)));
            axis image;
            colormap(plot_layout_cmap);
            title('Source Mask');            
        end

        % plot the sensor mask
        if record.use_sensor
            subplot(2, 2, 2), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(sensor.mask(x1:x2, y1:y2)));
            axis image;
            colormap(plot_layout_cmap);
            title('Sensor Mask');
        end

        % plot the sound speed distribution
        plot_zlim = double([min(c(:)), max(c(:))]);
        plot_rng = plot_zlim(2) - plot_zlim(1);
        if numel(c) == 1
            subplot(2, 2, 3), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(c)*ones(size(p(x1:x2, y1:y2))));
        else
            subplot(2, 2, 3), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(c(x1:x2, y1:y2)), [plot_zlim(1) - plot_spc*plot_rng, plot_zlim(2) + plot_spc*plot_rng]);
        end
        axis image;
        colormap(plot_layout_cmap);
        title('Sound Speed');

        % plot the density distribution
        plot_zlim = double([min(rho0(:)), max(rho0(:))]);
        plot_rng = plot_zlim(2) - plot_zlim(1);
        if numel(rho0) == 1
            subplot(2, 2, 4), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(rho0)*ones(size(p(x1:x2, y1:y2))));
        else
            subplot(2, 2, 4), imagesc(kgrid.y_vec(y1:y2)*scale, kgrid.x_vec(x1:x2)*scale, double(rho0(x1:x2, y1:y2)), [plot_zlim(1) - plot_spc*plot_rng, plot_zlim(2) + plot_spc*plot_rng]);
        end
        axis image;
        colormap(plot_layout_cmap);
        title('Density');  
        
        % add axis label
        xlabel(['(All axes in ' prefix 'm)']);        
    
    case 3
        
        % set the spacing as a fraction of the plot range
        plot_spc = 0.6;            
        
        % plot the initial pressure distribution
        if isfield(source, 'p0')
            figure;
            plot_zlim = double([min(source.p0(:)), max(source.p0(:))]);
            planeplot(scale*kgrid.x_vec(x1:x2), scale*kgrid.y_vec(y1:y2), scale*kgrid.z_vec(z1:z2), double(source.p0(x1:x2, y1:y2, z1:z2)), 'Initial Pressure: ', plot_zlim, prefix, plot_layout_cmap);
        end

        % plot c if heterogeneous
        if numDim(c) == 3
            figure;
            plot_zlim = double([min(c(:)), max(c(:))]);
            plot_rng = plot_zlim(2) - plot_zlim(1);            
            planeplot(scale*kgrid.x_vec(x1:x2), scale*kgrid.y_vec(y1:y2), scale*kgrid.z_vec(z1:z2), double(c(x1:x2, y1:y2, z1:z2)), 'c: ', [plot_zlim(1) - plot_spc*plot_rng, plot_zlim(2) + plot_spc*plot_rng], prefix, plot_layout_cmap);
        end

        % plot rho0 if heterogeneous
        if numDim(rho0) == 3    
            figure;
            plot_zlim = double([min(rho0(:)), max(rho0(:))]);
            plot_rng = plot_zlim(2) - plot_zlim(1);
            planeplot(scale*kgrid.x_vec(x1:x2), scale*kgrid.y_vec(y1:y2), scale*kgrid.z_vec(z1:z2),  double(rho0(x1:x2, y1:y2, z1:z2)), 'rho0: ', [plot_zlim(1) - plot_spc*plot_rng, plot_zlim(2) + plot_spc*plot_rng], prefix, plot_layout_cmap);
        end
        
end

% clean up unused variables
clear plot_ylim plot_rng plot_spc plot_layout_cmap