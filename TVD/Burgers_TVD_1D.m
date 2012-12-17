%% 1D Conservation Laws Solver using TVD
% Solving Invicid Burgers Equation
%
% u_t + f(u) = 0
%
% Where: f(u) = u^2/2
% 
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.12.17

clear all; close all; clc;

%% Parameters
     dx = 0.01;     % Spatial step size
    cfl = 0.8;      % Courant Number
 tStart = 0;        % Start time
   tEnd = 1;      % End time
IC_case = 3;        % {1} Gaussian, {2} Slope, {3} Triangle
limiter = 1;        % Options: 1(Vl), 2(Sb), 3(Mm), 4(koren)

%% Define our Flux function
    f = @(w) w.^2/2;

%% Discretization of Domain
      x = 1:dx:2;     % x grid
      nx = length(x);  % number of points

%% Initial Condition
% NOTE: In IC function, we can control the initial max and min value of u. 
     u0 = IC_iBurgers(x,3); 

% Load initial time
    time = tStart;
     
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
Fr = zeros(1,nx);       % Right flux    
Fl = zeros(1,nx);       % Left flux
%u_mp = zeros(1:nx-1)   % u values at the cell boundaries

while time < tEnd
 
% Update time step
    dt = cfl*dx/abs(max(u));    % time step size
  dtdx = dt/dx;                 % precomputed to save some flops
  time = time + dt;             % iteration actual time.

    % Compute fluxes at cell boundaries (middle points x_i+1/2)
    % Computing Flux to the Right: f(u_{i+1/2}) = 1/2(f(u_i)-f(u_i+1))
        
        u_mp = 1/2*(u(2:nx) - u(1:nx-1)); 
    % Compute left and right flux for every i-cell
    for i = 2:nx-1 % cells with two boundaries
        if u(i) >= 0;
            Fr(i) = f(u_mp(i));
        else
            Fr(i) = 0;
        end
        if u(i) < 0;
            Fl(i) = f(u_mp(i-1));
        else
            Fl(i) = 0;
        end
    end
    
    % Compute solution of next time step using Upwind
    u_next = u - dtdx * (Fr - Fl);

    % Periodic BC
    u_next(1)  = u(1);   % left boundary condition (fixed value)
    u_next(nx) = u(nx);  % right boundary condition (fixed value)
    
    % Update information
    u = u_next;
    
% Counter
it_count = it_count + 1;        % Update counter
  
end