%% 1D Boltzmann Equation
% Numerical solution of the boltzmann equation to recover Euler macroscopic
% continuum solution. By Manuel Diaz 2012.10.6
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
clc;  clear all;  close all;

%% Simulation Parameters
CFL    = 0.6;    % CFL condition
r_time = 0.01;   % Relaxation time
tEnd   = 0.2;    % End time
theta  = 0;      % for FD = -1, MB = 0, BE = +1
quad   = 1;      % for NC = 1 , GH = 2
method = 1;      % for TVD = 1, WENO3 = 2, WENO5 = 3
input  = 2;      % Reimann IC case

%% Space Discretization
nx  = 40;                       % Desided number of points in our domain
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-10,10];  % range: a to b
    nv = 40;        % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    nv = 20;        % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
     k = 1;         % quadrature constant.
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual nv value will be computed using 'lenght' vector function:
nv = length(v); 

% Initialize Arrays
% u = zeros(nv,nx);     t = zeros(nv,nx);   r = zeros(nv,nx); 
n = zeros(nv,nx);     p = zeros(nv,nx);   
v = repmat(v,1,nx);   w = repmat(w,1,nx);


%% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
    [r0,u0,t0] = reimann_IC1d(x,input);
    
% Using Discrete Ordinate Method:
    r = repmat(r0,nv,1); u = repmat(u0,nv,1); t = repmat(t0,nv,1);

% Compute distribution IC of our mesoscopic method by assuming the equilibrium 
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
    f0 = f_equilibrium_1d(r,u,v,t,theta);

% Plot initial Equilibrium function:
    surface(f0); grid on;
    xlabel('x - Spatial Domain'); ylabel('v - Velocity Space');
    zlabel('f - Probability');

% Compute Initial Macroscopic Momemts:
    %[n,j_x,E] = macromoments1d(k,w,f,v);
    
%% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evonve the same problem for
% every mesoscopic velocity.
dt = dx*CFL/max(v(:,1)); 
dtdx = dt/dx;  % precomputed to save someflops

% time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the classic
% transport Boltzmann equation will resemble to a pure advection equation.
% Thus WENO, TVD, DG or CPR can be used easyly to compute evolution of the
% information inside the domain: 

switch method
    
    case{1} % TVD 0(h^2)
        % Using discrete ordinate method (discrete and constant velocity
        % values in phase-space domain)
        a = v(:,1);
        
        % Load initial condition
        f = f0;
        
        % Initialize vector variables
        f_next = zeros(1,nx);
        
        for tsteps = time
            % for visualization
            figure(1)
            h = plot(x,r(1,:)); axis([0,1,0,0.5])
        
            % compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(r,u,v,t,theta);
                  
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % Compute the smoothness factors, r(j), from data, u(j).
                [r] = theta1d(f(i,:),a(i));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(f(i,:),a(i),dtdx,phi);

                % Compute next time step
                f_next(i,:) = f(i,:) - dtdx*(F_right - F_left) ...
                    + (dt/r_time)*(f_eq(i,:)-f(i,:));

                % BC
                f_next(i,1) = f_next(i,2);
                f_next(i,nx) = f_next(i,nx-1);

                % UPDATE info
                f(i,:) = f_next(i,:);
            end
            
            % Compute macroscopic moments
            [n,j_x,E] = macromoments1d(k,w,f,v);
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing chalenge)
            [r,u,t,p] = macroproperties1d(n,j_x,E,nx,nv,theta);
            
            %fprintf('%2.4f\n',tsteps)
            drawnow
        end
        
    case{2} % WENO k = 3 i.e. O(h^5) 
        
    otherwise
        error('Order must be between 1 and 2');
end

%% Plot results
%plot(x,u);