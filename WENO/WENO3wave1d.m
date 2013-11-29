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
% Notes: Basic Scheme Implementation without RK integration method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

%% Parameters
    cfl = 0.20;  % Courant Number
     nx = 80;    % number of cells
   tEnd = 2*pi;  % End time
	 ic = 1;     % {1} Gaussian, {2} SquareJump, {3} Riemann, {4} Sine 
BC_type = 1;     % {1} Dirichlet, {2} Neumann, {3} Periodic
fluxsplit = 2;   % {1} Godunov, {2} Global LF, {3} Local LF


%% Define our Flux function
fluxtype = 'linear';

switch fluxtype
    case 'linear'
        c = 1; flux = @(w) c*w; % scalar advection speed
        dflux = @(w) c*ones(size(w));
    case 'nonlinear'
        flux = @(w) w.^2/2;
        dflux = @(w) w;
end

%% Discretization of Domain
% Depending on the size of the degree of the WENO stencil we are using then
% 2 or more ghost cells have to be introduced
a=0; b=2*pi; x=linspace(a,b,nx); dx=(b-a)/nx; 

%% Initial Condition
u0 = IC(x,ic); 
          
%% Solver Loop
% Beause the max slope, f'(u) = u, may change as the time steps progress
% we cannot fix the size of the time steps. Therefore we need to recompute
% the size the next time step, dt, at the beginning of each iteration to
% ensure stability.

% WENO Interpolation arrays
ur = zeros(1,nx);	% u values at x_{i+1/2}^{-} (right)
ul = zeros(1,nx);	% u values at x_{i-1/2}^{+} (left)

% Build flux maps 
mapr = [1,1:nx];
mapl = [1:nx,1];

% load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    % Update
    dt = cfl*dx/abs(max(u));  % time step
    dtdx = dt/dx;             % precomputed
    t = t+dt;	% actual time.
    it = it+1;  % iteration counter
    
    % Plot solution
    plotrange = 3:nx-2;
    plot(x,u0,'-x',x(plotrange),u(plotrange),'.'); 
    axis([a,b,min(u0)-0.1,max(u0)+0.1])
    xlabel('X-Coordinate [-]'); ylabel('U-state [-]'); 
    title 'Invicid Burgers using WENO scheme';
    
    % Left and Right weno reconstructions (this can be done in parallel)
    for i = 3:nx-2
        sr = i-2:i+2;   % WENO stencil range
        [ur(i),ul(i)] = WENO3_1d_flux(u(sr),u(sr));    
    end
    %plot(x,ur(mapr),'o',x,ul(mapl),'+',x,u,'-.'); % ok!
    
    % Flux Spliting
    alpha = max(dflux(u));
    h = 0.5*( flux(ul(mapl)) + flux(ur(mapr))... % the average 
        -alpha*( ur(mapr) - ul(mapl) )); % - the jump
      
    % Update Next time step
    u_next = u - dtdx *(h(2:end)-h(1:end-1));
    
    % Apply BCs
    u_next = WENO3_1d_BCs(u_next,BC_type,nx);
    
    % Update Info
    u = u_next;
    
    % Update plot
    drawnow
end