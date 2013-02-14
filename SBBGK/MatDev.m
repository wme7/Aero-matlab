function [DuDt] = MatDev(u,u_old,c,dtdx)
% This function evaluates the material Derivate:  Du/Dt = du/dt + c du/dx
% in this function, u and u_next must be vector arrays of the same size.

strategy = 2;
switch strategy
    case{1} % Upwind
        [cu_left,cu_right] = Upwindflux1d(u,c); % Compute Upwind Fluxes
        DuDt = (u - u_old) + dtdx.*(cu_right-cu_left);
    case{2} % TVD
        [r] = theta1d(u,c); % Compute the smoothness factors, r(j), from data, u(j).
        [phi] = fluxlimiter1d(r,1); % Compute the Flux Limiter, using Van Leer Limiter
        [cu_left,cu_right] = TVDflux1d(u,c,dtdx,phi); % Compute TVD Fluxes
        DuDt = (u - u_old) + dtdx.*(cu_right-cu_left);
end