function [DuDt] = MatDev(u,u_old,c,dtdx)
% This function evaluates the material Derivate:  Du/Dt = du/dt + c du/dx
% in this function, u and u_next must be vector arrays of the same size.

[cu_left,cu_right] = Upwindflux1d(u,c);
DuDt = (u - u_old) + dtdx.*(cu_right-cu_left);