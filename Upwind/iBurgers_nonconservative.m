%% Non-conservative method for Non-linear eq. in 1D 
% to solve the scalar advection equation:
%
% du/dt + u du/dx = 0 ; Burgers' equation
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2013.10.12

clear all; close all; clc;

%% Parameters
     dx = 0.01; % Spatial step size
    cfl = 0.35;  % Courant Number
  t_end = 0.2;    % End time

%% Discretization of Domain
x = 1:dx:2;     % x grid
n = length(x);  % number of points

%% Initial Condition
%u0 = IC_iBurgers(x,2);
u0 = ones(size(x));
xc = 1.5; rhs = find(x<xc);
u0(rhs) = 2;

%% Double-Sided Upwind

% Load IC
u_upwind = u0;

% Main loop
t = 0;
while t < t_end
    % update time step
    dt = cfl*dx/max(abs(u_upwind)); % time step size
    dtdx = dt/dx;	% precomputed to save some flops
    t = t + dt;
    
    % Plot solution
    plot(x,u_upwind,'.')
    
    for j = 2:n-1
    % Define positive and negative wave velocites
    a_p = max(0,u_upwind(j)); % a^{+}
    a_m = min(0,u_upwind(j)); % a^{-}
    % Next time step
    u_next(j) = u_upwind(j) ...
        - a_p*dt/dx*(u_upwind(j) - u_upwind(j-1)) ...
        - a_m*dt/dx*(u_upwind(j+1) - u_upwind(j));
    end
    u_upwind = u_next;
    % BC
    u_upwind(1) = u_next(2);
    u_upwind(n) = u_next(n-1);
    
    % Update plot
    drawnow;
end
u_next = zeros(1,n);
