function [h_left,h_right] = WENO5_1d_driver(vp,vn)
% WENO3 Reconstruction Driver

% Number of cells
nx = length(vp);

% Allocate variables
hn = zeros(size(vn));    
hp = zeros(size(vp));

% Reconstruc h+ and h- fluxes
for i = 3:nx-2
    xr = i-2:i+2;   % x-range of cells
    [hn(i),hp(i-1)] = WENO5_1d_flux(vp(xr),vn(xr));    
end
h = sum([hn;hp]); 

% Formulate Left and Right fluxes, equiv: %h(i) - h(i-1)
h_right = [h(1:end-2),0,0];
h_left = [0,0,h(2:end-1)];