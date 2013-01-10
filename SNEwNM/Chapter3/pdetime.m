function ft=pdetime(u)
% PDETIME 
% FT = PDETIME(U) is the nonlinear residual for time-dependent 
% problem in Chapter 3.
% This code has the zero boundary conditions, the right
% side, and the preconditioner built in.
%
% The equation is u_t = (u_xx + u_yy) - 20*u*(u_x + u_y)  + f
%
% Where f is the right hand side of the steady-state example
%
%
global rhsf uold dt
%
% Left preconditioned nonlinear residual for implicit Euler discretization.
%
v=20*u.*(dxmf(u)+dymf(u))-rhsf;
%
%
ft=dt*u + fish2d(u - uold + dt*v);

