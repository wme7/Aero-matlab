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
     nx = 100;  % Spatial step size
    cfl = 0.1;	% Courant Number
   tEnd = 0.2;	% End time

%% Define our Flux function
    c = 1; flux = @(w) c.*w;
% and the Derivate of the flux function
dflux = @(w) a.*ones(size(w));
      
%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 2 or more ghost cells have to be introduced
a=0; b=2*pi; x=linspace(a,b,nx); dx=(b-a)/nx; 

%% Initial Condition
u0 = IC(x,6); % IC: {1}Gaussian, {2}Sine, {3}CSine, {4}Riemann, {5}Tanh, {6}SquareJump
    
%% Main Loop
% Load initial condition into the domain
    u = u0;  

% Time discretization
	dt = cfl*dx/abs(c); % time step size
  dtdx = dt/dx;	% precomputed to save some flops
     t = 0:dt:tEnd;

% WENO interpolations vector
 u_next = zeros(1,nx);
     ur = zeros(1,nx-1);	% u values at the right cell boundaries
     ul = zeros(1,nx-1);	% u values at the right cell boundaries

tic
for kk = t
    if c > 0;
        % Upwinding WENO flux
        for i = 4:nx-3
            su = u(i-3:i+3);
            ur(i) = WENO7_upwind(su);
        end

        % Compute solution of next time step using WENO Upwind
        for i = 4:nx-3
            u_next(i) = u(i) - c*dtdx*(ur(i) - ur(i-1));
        end

    else % c < 0;
        % Downwinding WENO flux
        for i = 4:nx-3
            su = u(i-3:i+3);
            ul(i) = WENO7_downwind(su);
        end

        % Compute solution of next time step using WENO Upwind
        for i = 4:nx-3
            u_next(i) = u(i) - c*dtdx*(ul(i+1) - ul(i));
        end
    end

    % Update BCs: All Neumann BC's
   	u_next = WENO5_1d_BCs(u_next,2,nx); % 2: Neumann BC
    
    % UPDATE info
    u = u_next;
end
toc

%% Exact Solution
% Displacement:
x_add = floor(c*tEnd/dx);
% Computing exact solution:
u_exact = ones(1,nx);
 x_1 = ceil(2*nx/5) + x_add;
 x_2 = ceil(3*nx/5) + x_add;
u_exact(x_1:x_2) = 2;

%% Plot Results
hold on
plot(x,u,'.'); plot(x,u_exact,'k'); 
xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); 
ylim([0.5,2.5]); title 'WENO5';
hold off