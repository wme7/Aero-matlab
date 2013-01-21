%% Weighted Essentially non-Oscilaroty (WENO) Routine
% To solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
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
     dx = 0.01;  % Spatial step size
    cfl = 0.80;  % Courant Number
 tStart = 0.00;  % Start time
   tEnd = 0.15;  % End time
IC_case = 4;     % {1} Gaussian, {2} Slope, {3} Triangle, {4} Sine
flxtype = 3;     % {1} Roe, {2} LF, {3} LLF, {4} Upwind <-non-conservative!

%% Define our Flux function
     f = @(w) w.^2/2;
% and the Derivate of the flux function
    df = @(w) w;

%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 1 or more ghost cells have to be introduced
    weno = 3;   gcells = weno-1;
    a = 0 - gcells*dx;
    b = 1 + gcells*dx;
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
h = zeros(1,nx);        % Flux values at the cell boundaries
hn = zeros(1,nx);       % Flux values at u_i+1/2 (-)
hp = zeros(1,nx);       % Flux values at u_i+1/2 (+)

while time <= tEnd
    % Plot Evolution
    range = 3:nx-2;
    plot(x(range),u(range),'.'); 
    axis([x(1) x(end) min(u0)-0.1 max(u0)+0.1])
    xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); 
    title 'Invicid Burgers using WENO scheme';
    
    % Update time step
    dt   = cfl*dx/abs(max(u));  % time step size
    dtdx = dt/dx;               % precomputed to save some flops
    time = time + dt;           % iteration actual time.
    
    % Numerical Global Lax Friedrichs Flux Spliting
    %alpha = max(abs(df(u)));
    %fp = 0.5*( f(u) + alpha*u );
    %fn = 0.5*( f(u) - alpha*u );
    
    % Compute WENO Flux values at cells interfaces
    for i = 3:nx-2
        [hn(i),hp(i)] = WENO3_1D_flux(f(u(i-2:i+2)));
    end
    h = hn + hp;
    
    % Numerical Flux Spliting
    
    
    % Compute solution of next time step using WENO Upwind
    for i = 4:nx-3
    u_next(i) = u(i) - dtdx * (h(i) - h(i-1));
    end

    % Periodic BC
%     u_next(nx-2) = u_next(4);    % right boundary condition (periodic)
%     u_next(nx-1) = u_next(5);    % right boundary condition (periodic)
%     u_next( nx ) = u_next(6);    % right boundary condition (periodic)
%     u_next(1)    = u_next(nx-5); % left boundary condition (periodic)
%     u_next(2)    = u_next(nx-4); % left boundary condition (periodic)
%     u_next(3)    = u_next(nx-3); % left boundary condition (periodic)
        
    % Dirichlet BC
    fix = 0.5;
    u_next(1)    = fix; % left boundary condition (fixed value)
    u_next(2)    = fix; % left boundary condition (fixed value)
    u_next(3)    = fix; % left boundary condition (fixed value)
    u_next(nx-2) = fix; % right boundary condition (fixed value)
    u_next(nx-1) = fix; % right boundary condition (fixed value)
    u_next( nx ) = fix; % right boundary condition (fixed value)
    
    % UPDATE Info
    u = u_next;
    
    % Counter
    it_count = it_count + 1; % Update counter
    
    % draw iteration
    drawnow
end