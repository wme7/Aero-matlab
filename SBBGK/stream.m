function [u_next] = stream(u,c,dtdx,strategy)
% This function solver the following 1D nonlinear conservation equation,
% $du/dt + df(u)/dx = 0$ 
% coded by Manuel Diaz, NTU, 2013.03.29.

% Advection of Information
switch strategy
    case{1} % Upwind, O(h)
        [f_left,f_right] = Upwindflux1d(u,c); % Compute Upwind Fluxes
        u_next = u - dtdx.*(f_right-f_left);
        
        % BCs: Neumann
        u_next(1) = u_next(2);      % left boundary condition
        u_next(end) = u_next(end-1);  % right boundary condition
        
    case{2} % TVD, 0(h^2)
        [r] = theta1d(u,c); % Compute the smoothness factors, r(j)
        [phi] = fluxlimiter1d(r,1); % Flux Limiter: {1} Van Leer Limiter
        [f_left,f_right] = TVDflux1d(u,c,dtdx,phi); % Compute TVD Fluxes
        u_next = u - dtdx.*(f_right-f_left);
        
        % BCs: Neumann
        u_next(1) = u_next(2);      % left boundary condition
        u_next(end) = u_next(end-1);  % right boundary condition
        
    case{3} % WENO3, 0(h^5)
        [vp,vn] = WENO_scalarfluxsplit(c.*u); % Slip Positive and Negative Scalar Fluxes
        [f_left,f_right] = WENO3_1d_driver(vp,vn); % Compute WENO3 Fluxes
        u_next = u - dtdx.*(f_right-f_left);
        
        % BCs: Neumann
        u_next(end-2) = u_next(end-3); % right boundary condition
        u_next(end-1) = u_next(end-3); % right boundary condition
        u_next( end ) = u_next(end-3); % right boundary condition
        u_next(3)    = u_next(4);    % left boundary condition
        u_next(2)    = u_next(4);    % left boundary condition
        u_next(1)    = u_next(4);    % left boundary condition
        
    case{4} % WENO5, 0(h^9)
        [vp,vn] = WENO_scalarfluxsplit(c.*u); % Slip Positive and Negative Scalar Fluxes
        [f_left,f_right] = WENO5_1d_driver(vp,vn); % Compute WENO3 Fluxes
        u_next = u - dtdx.*(f_right-f_left);
        
        % BCs: Neumann
        u_next(end-3) = u_next(end-4); % right boundary condition
        u_next(end-2) = u_next(end-4); % right boundary condition
        u_next(end-1) = u_next(end-4); % right boundary condition
        u_next( end ) = u_next(end-4); % right boundary condition
        u_next(1)    = u_next(5);    % left boundary condition
        u_next(2)    = u_next(5);    % left boundary condition
        u_next(3)    = u_next(5);    % left boundary condition
        u_next(4)    = u_next(5);    % left boundary condition
        
    otherwise
        error('stream case not yet available');
        
end

