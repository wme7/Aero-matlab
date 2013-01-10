function ft=ftime(u)
% FTIME 
% Nonlinear residual for time-dependent problem in Chapter 2.
% This code has the zero boundary conditions built in.
% The time step and solution are passed as globals.
%
global uold dt
%
% d2u is the numerical negative second derivative.
%
n=length(u); h=1/(n+1);
d2u=2*u;
d2u(1:n-1)=d2u(1:n-1)-u(2:n);
d2u(2:n)=d2u(2:n)-u(1:n-1);
d2u=d2u/(h^2);
%
% Nonlinear residual for implicit Euler discretization.
%
ft=(u - uold) - dt * (exp(u) - d2u);

