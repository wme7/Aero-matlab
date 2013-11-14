%% 1D Conservation Laws Solver using Upwind
% To solve Invicid Burgers Equation
%
%   $u_t + f(u) = 0$ for x in [a,b]
%
% where: $f(u) = u^2/2$
% 
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.12.17

clear all; clc; close all;

%% Parameters
     dx = 1/200; % Spatial step size
    cfl = 1/2;   % Courant Number
   tEnd = 2.00;  % End time
IC_case = 4;     % {1} Gaussian, {2} Slope, {3} Triangle, {4} Sine
limiter = 1;     % Options: 1(Vl), 2(Sb), 3(Mm), 4(koren)
flxtype = 4;     % {1} Godunov, {2} Roe, {3} LF, {4} LLF, {5} Upwind {not conservative}

%% Define our Flux function
     f = @(w) w.^2/2;
% and the Derivate of the flux function
    df = @(w) w;

%% Discretization of Domain
      x = 0:dx:2;     % x grid
      nx = length(x); % number of points

%% Initial Condition
% NOTE: In IC function, we can control the initial max and min value of u. 
     u0 = IC(x,IC_case); 

% Load initial time
    t = 0;
     
% Load Initial Condition
    u = u0;
    
%% Main Loop
% Beause the max slope, f'(u) = u, is changing as the time steps progress
% we cannot fix the size of the time steps. Therefore we need to recompute
% the size the next time step, dt, at the beginning of each iteration to
% ensure stability in our routine. 

% Initialize vector variables 
it_count = 0;           % Iteration counter
u_next = zeros(1,nx);   % u in next time step
h = zeros(1,nx);        % Flux values at the cell boundaries

while t < tEnd
    % Plot Evolution
    plot(x,u); axis([x(1) x(end) min(u0)-0.1 max(u0)+0.1])
    
    % Update time step
    dt   = cfl*dx/abs(max(u));  % time step size
    dtdx = dt/dx;               % precomputed to save some flops
    t = t + dt;           % iteration actual time.

    % Compute fluxes at cell boundaries (middle points x_i+1/2)
    h = flux1d(f,df,u,flxtype);
    
    % Compute solution of next time step using Upwind
    for i = 2:nx-1
        u_next(i) = u(i) - dtdx * (h(i) - h(i-1));
    end

    % Periodic BC
    u_next(1)  = u(nx);   % left boundary condition (periodic)
    u_next(2)  = u(1);    % left boundary condition (periodic)
    u_next(nx) = u(nx-1); % right boundary condition (fixed value)
    
    % Dirichlet BC
    %u_next(1)  = u(1);    % left boundary condition (fixed value)
    %u_next(nx) = u(nx);   % right boundary condition (fixed value)
    
    % Update information
    u = u_next;
    
% Counter
it_count = it_count + 1; % Update counter

% draw iteration
drawnow

end
