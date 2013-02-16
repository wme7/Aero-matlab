function [vp,vn] = WENO_scalarfluxsplit(v)
% WENO flux spliting subroutine.
% OUTPUT:
% * vp: positive flux, v^{+}, which is intended for reconstructing f_{i+1/2}^{-}
% * vn: negative flux  v^{-}, which in intended for reconstructing f_{i+1/2}^{+}

% Godunov Flux Spliting
vp = 0.5*(v + abs(v)); %flux^{+}
vn = 0.5*(v - abs(v)); %flux^{-}
