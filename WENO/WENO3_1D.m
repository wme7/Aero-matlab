%% Weighted Essentially non-Oscilaroty (WENO) Routine
% to solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
%
% where $F(u) = a*u$
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.20

clear all; close all; clc;

%% Parameters
      a = -1.10;     % Scalar velocity in x direction
     dx = 0.01;     % Spatial step size
    cfl = 0.62;     % Courant Number
     dt = cfl*dx/abs(a); % time step size
   dtdx = dt/dx;    % precomputed to save some flops
   tEnd = 0.25;     % End time

%% Define our Flux function
     f = @(w) a.*w;
% and the Derivate of the flux function
    df = @(w) a.*ones(size(w));
      
%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 1 or more ghost cells have to be introduced
    r = 3;              % for WENO r=3
    gcells = r - 1;     % gosht cells to add
    xa = 1 - gcells*dx; % a
    xb = 2 + gcells*dx; % b
    x = xa:dx:xb;       % x grid
    nx = length(x);      % number of points
    
%% Time Discretization
t = 0:dt:tEnd;

%% Initial Condition
u_0 = ones(1,nx);
 x_1 = ceil(2*nx/5);
 x_2 = ceil(3*nx/5);
u_0(x_1:x_2) = 2;

%% Main Loop
% Load initial condition into the domain
u = u_0; 

% Initialize vector variables 
u_next = zeros(1,nx);
h = zeros(1,nx-1);      % Flux values at the cell boundaries

tic
for kk = t
    if a > 0;
        % Upwinding WENO flux
        for i = 3:nx-2
            fu = f(u(i-2:i+2));
            h(i) = WENO_upwind(fu);
        end

        % Compute solution of next time step using WENO Upwind
        for i = 4:nx-3
            u_next(i) = u(i) - dtdx * (h(i) - h(i-1));
        end

    else % a < 0;
        % Downwinding WENO flux
        for i = 3:nx-2
            fu = f(u(i-2:i+2));
            h(i) = WENO_downwind(fu);
        end

        % Compute solution of next time step using WENO Upwind
        for i = 4:nx-3
            u_next(i) = u(i) - dtdx * (h(i+1) - h(i));
        end
    end

    % Update BCs: All Neumann BC's
   	u_next = WENO3_1d_BCs(u_next,2,nx); % 2: Neumann BC
    
    % UPDATE info
    u = u_next;
end
toc

%% Exact Solution
% Displacement:
x_add = floor(a*tEnd/dx);
% Computing exact solution:
u_exact = ones(1,nx);
 x_1 = ceil(2*nx/5) + x_add;
 x_2 = ceil(3*nx/5) + x_add;
u_exact(x_1:x_2) = 2;

%% Plot Results
hold on
plot(x,u,'.'); plot(x,u_exact,'k'); xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); ylim([0.5,2.5]); title 'WENO3';
hold off