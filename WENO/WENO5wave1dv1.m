%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Solving 1-D wave equation with WENO
%                Weighted Essentially non-Oscilaroty
%
%                du/dt + df/dx = 0,  for x \in [a,b]
%                  where f = f(u): linear/nonlinear
%
%             codedby Manuel Diaz, manuel.ade'at'gmail.com 
%              Institute of Applied Mechanics, 2012.08.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: Jiang & Shu JCP. vol 126, 202-228 (1996)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Finite Difference Scheme without RK integration method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

%% Parameters
      cfl = 0.30;   % Courant Number
       nx = 80;     % number of cells
     tEnd = 2*pi;   % End time
  BC_type = 2;      % {1} Dirichlet, {2} Neumann, {3} Periodic
fluxsplit = 2;      % {1} Godunov, {2} Global LF, {3} Local LF

%% Define our Flux function
fluxtype = 'linear';

switch fluxtype
    case 'linear'
        i = 1; flux = @(w) i*w; % scalar advection speed
        dflux = @(w) i*ones(size(w));
    case 'nonlinear'
        flux = @(w) w.^2/2;
        dflux = @(w) w;
end

%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 2 or more ghost cells have to be introduced
a=0; b=2*pi; x=linspace(a,b,nx); dx=(b-a)/nx; 

%% Initial Condition
u0 = IC(x,5); % IC: {1}Gaussian, {2}Sine, {3}CSine, {4}Riemann, {5}Tanh
          
%% Solver Loop
% Beause the max slope, f'(u) = u, may change as the time steps progress
% we cannot fix the size of the time steps. Therefore we need to recompute
% the size the next time step, dt, at the beginning of each iteration to
% ensure stability.

% load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    % Update
    dt = cfl*dx/abs(max(u));  % time step
    dtdx = dt/dx;             % precomputed to save some flops
    t = t+dt;	% actual time.
    it = it+1;  % iteration counter
    
    % Plot solution
    plotrange = 3:nx-2;
    plot(x,u0,'-x',x(plotrange),u(plotrange),'.'); 
    axis([a,b,min(u0)-0.1,max(u0)+0.1])
    xlabel('X [-]'); ylabel('U [-]'); title 'Invicid Burgers using WENO';
    
    % Flux Spliting
    [vp,vn] = fluxSplit(u,flux,dflux,fluxsplit);

    % Reconstruct Fluxes values at cells interfaces
    [h_left,h_right] = WENO5_1d_driver(vp,vn);
    
    % Update Next time step
    u_next = u - dtdx *(h_right - h_left);

    % Apply BCs
    u_next = WENO5_1d_BCs(u_next,BC_type,nx);
        
    % Update Info
    u = u_next;
    
    % Update plot
    drawnow
end