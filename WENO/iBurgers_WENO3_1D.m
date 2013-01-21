%% Weighted Essentially non-Oscilaroty (WENO) Routine
% To solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0  for x in [a,b]
%
% Where can be: 
%   $F(u) = a*u$
%   $F(u) = u.^2/2$
%
% Ref:
% Jiang & Shu JCP. vol 126, 202-228 (1996)
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.20

clear all; clc; close all;

%% Parameters
     dx = 0.02;  % Spatial step size
    cfl = 0.80;  % Courant Number
 tStart = 0.00;  % Start time
   tEnd = 3.50;  % End time
IC_case = 5;     % {1} Gaussian, {2} Slope, {3} Triangle, {4} Sine {5} Riemann
bc_type = 1;     % {1} Dirichlet, {2} Neumann, {3} Periodic

%% Define our Flux function
     f = @(w) w.^2/2;
% and the Derivate of the flux function
    df = @(w) w;

%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 1 or more ghost cells have to be introduced
    r = 3;              % for WENO r=3
    gcells = r-1;       % gosht cells to add
    a = 1 - gcells*dx;  % a
    b = 2 + gcells*dx;  % b
    x = a:dx:b;         % x grid
    nx = length(x);     % number of points

%% Initial Condition
% NOTE: In IC function, we can control the initial max and min value of u. 
     u0 = IC_iBurgers(x,IC_case); 

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
h = zeros(2,nx-1);      % Flux values at the cell boundaries
fn = zeros(1,nx-1);     % Flux values at u_i+1/2 (-)
fp = zeros(1,nx-1);     % Flux values at u_i-1/2 (+)

while time <= tEnd
    % Plot Evolution
    range = 3:nx-2;     % plot and solution range
    plot(x(range),u(range),'.'); 
    axis([x(1) x(end) min(u0)-0.1 max(u0)+0.1])
    xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); 
    title 'Invicid Burgers using WENO scheme';
    
    % Update time step
    dt   = cfl*dx/abs(max(u));  % time step size
    dtdx = dt/dx;               % precomputed to save some flops
    time = time + dt;           % iteration actual time.
    
    % Compute WENO Flux values at cells interfaces
    for i = 3:nx-2
        [fn(i),fp(i+1)] = WENO3_1D_flux(f(u(i-2:i+2)));
    end
    fluxes = [fn;fp];
    
    % Numerical Global Lax Friedrichs Flux Spliting
    du = (u(2:nx)-u(1:nx-1));  alpha = max(abs(df(u))); 
    h = 0.5 * (fluxes(2,:) + fluxes(1,:) - alpha * du);
        
    % Compute solution of next time step using WENO Upwind
    for i = 4:nx-3
        u_next(i) = u(i) - dtdx * (h(i) - h(i-1));
    end

    % WENO BC's
    u_next = WENO3_1D_BCs(u_next,bc_type,nx);
    
    % UPDATE Info
    u = u_next;
    
    % Counter
    it_count = it_count + 1; % Update counter
    
    % draw iteration
    drawnow
end