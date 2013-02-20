function [h_left,h_right] = WENO5_1d_driver(vp,vn)
% WENO3 Reconstruction Driver

% Number of cells
nx = length(vp);

% Allocate variables
hn = zeros(size(vn));    
hp = zeros(size(vp));

% Reconstruc h+ and h- fluxes
for i = 4:nx-3
    xr = i-3:i+3;   % x-range of cells
    [hn(i),hp(i-1)] = WENO5_1d_flux(vp(xr),vn(xr));    
end
h = sum([hn;hp]); 

% Formulate Left and Right fluxes, equiv: %h(i) - h(i-1)
h_right = [h(1:end-3),0,0,0];
h_left = [0,0,0,h(3:end-1)];