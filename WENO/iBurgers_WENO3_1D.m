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
    cfl = 0.20;  % Courant Number
     dx = 0.02;  % Spatial step size
 tStart = 0.00;  % Start time
   tEnd = 3.15;  % End time
IC_case = 1;     % {1} Gaussian, {2} Slope, {3} Triangle, {4} Sine {5} Riemann
bc_type = 3;     % {1} Dirichlet, {2} Neumann, {3} Periodic
fluxsplit = 2;   % {1} Godunov, {2} Global LF, {3} Local LF
fluxtype = 2;    % {1} Scalar, {2}, Burgers

%% Define our Flux function
switch fluxtype
    case{1}
        a = 0.5;   % scalar advection speed
        f = @(w) a*w;
        df = @(w) a*ones(size(w));
    case{2}
        f = @(w) w.^2/2;
        df = @(w) w;
    otherwise
        error('case not available')
end

%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 1 or more ghost cells have to be introduced
    k = 2;              % for WENO3 k=2
    gcells = k;         % gosht cells to add
    a = 1 - gcells*dx;  % a
    b = 3 + gcells*dx;  % b
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
%h = zeros(2,nx-1);      % Flux values at the cell boundaries
hn = zeros(1,nx-1);     % Flux values at x_i+1/2 (-)
hp = zeros(1,nx-1);     % Flux values at x_i-1/2 (+)

while time <= tEnd
    % Plot Evolution
    range = 3:nx-2;     % plot and solution range
    plot(x(range),u(range),'.'); 
    axis([x(1) x(end) min(u0)-0.1 max(u0)+0.1])
    xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); 
    title 'Invicid Burgers using WENO scheme';
    
    % Update time step
    dt   = cfl*dx/abs(max(u));  % size of time step
    dtdx = dt/dx;               % precomputed to save some flops
    time = time + dt;           % iteration actual time.
    
    % Flux Spliting
    [vp,vn] = WENO_fluxsplit(u,f,df,fluxsplit);

    % Reconstruct Fluxes values at cells interfaces
    for i = 3:nx-2
        xr = i-2:i+2; % x-range of cells
        %xr = i-3:i+3; % x-range of cells
        [hn(i),hp(i-1)] = WENO3_1d_flux(vp(xr),vn(xr));
        %[hn(i),hp(i-1)] = WENO5_1d_flux(vp(xr),vn(xr));
    end
    h = sum([hn;hp]);
                
    % Compute next time step
    for i = 4:nx-3
        u_next(i) = u(i) - dtdx * (h(i) - h(i-1));
    end

    % Apply BCs
    u_next = WENO3_1d_BCs(u_next,bc_type,nx);
    %u_next = WENO5_1d_BCs(u_next,bc_type,nx);
    
    % UPDATE Info
    u = u_next;
    
    % Counter
    it_count = it_count + 1; % Update counter
    
    % draw iteration
    drawnow
end