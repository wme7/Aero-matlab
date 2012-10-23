%% 2D Boltzmann Equation
% Numerical solution of the boltzmann equation to recover Euler macroscopic
% continuum solution. By Manuel Diaz 2012.10.06
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
clc;  clear all;  close all;

%% Simulation Parameters
CFL    = 0.6;   % CFL condition
r_time = 0.01;   % Relaxation time
tEnd   = 0.2;    % Out time
theta  = 0;      % for FD = -1, MB = 0, BE = +1
quad   = 2;      % for NC = 1 , GH = 2
method = 1;      % for TVD = 1, WENO3 = 2, WENO5 = 3
input  = 5;      % Reimann IC case

%% Space Discretization
nx = 20; ny = 20;
[X,dx,Y,dy] = grid2d(0,1,nx,0,1,ny);
[x,y] = meshgrid(X,Y);

%% Velocity Discretization:
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-10,10];  % range: a to b
    lv = 40;        % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),lv,5); % cotes Degree 5

    case{2} % Gauss Hermite Quadrature:
    lv = 20;        % nodes desired (the actual value)
    [v,w] = GaussHermite(lv); % for integrating range: -inf to inf
     k = 1;         % quadrature constant.
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual nv value will be computed using 'length' vector function:
nv = length(v); 

% Using Discrete Ordinate Method:
    v1 = v;  v1 = repmat(v1,[1,nv]); v1 = repmat(v1,[1,1,ny,nx]);
    v2 = v'; v2 = repmat(v2,[nv,1]); v2 = repmat(v2,[1,1,ny,nx]);
    w1 = v;  w1 = repmat(w1,[1,nv]); w1 = repmat(w1,[1,1,ny,nx]);
    w2 = v'; w2 = repmat(w2,[nv,1]); w2 = repmat(w2,[1,1,ny,nx]); 

% Initialize Arrays
    f  = zeros(nv,nv,ny,nx);    %n  = zeros(ny,nx);  
    u1 = zeros(nv,nv,ny,nx);    u2 = zeros(nv,nv,ny,nx);
    r  = zeros(nv,nv,ny,nx);    t  = zeros(nv,nv,ny,nx);    

%% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
    [r0,u1_0,u2_0,t0] = reimann_IC2d(x,y,input);
    
% Using Discrete Ordinate Method: (Replicating Data for easyer arrays Ops)
for j = 1:ny
    for i = 1:nx
    r(:,:,j,i)  = r0(j,i);      u1(:,:,j,i) = u1_0(j,i);
    u2(:,:,j,i) = u2_0(j,i);    t(:,:,j,i)  = t0(j,i);
    end
end

% Compute distribution IC of our mesoscopic method by assuming the equilibrium
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
%    f0 = f_equilibrium_2d(r,u1,u2,v1,v2,t,theta);

% Plot initial Equilibrium function:
%      surface(f0); grid on;
%      xlabel('x - Spatial Domain'); ylabel('v - Velocity Space');
%      zlabel('f - Probability');

% Compute Initial Macroscopic Momemts:
%    [n,j_x,j_y,E] = macromoments2d(k,w1,w2,f0,v1,v2);
    
%% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.
dt = min(dx,dy)*CFL/max(v); 
dtdx = dt/dx;  % precomputed to save someflops
dtdy = dt/dy;  % precomputed to save someflops

% time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the classic
% transport Boltzmann equation will resemble to a pure advection equation.
% Thus we use WENO, TVD, DG or CPR can be used to compute evolution of the
% information inside the domain: 

switch method
    
    case{1} % TVD 0(h^2)
        % Using Discrete Ordinate Method (Discrete and constant velocity
        % values in phase-space domain)
        a = v1(:,:,1,1);
        b = v2(:,:,1,1);
        
        % Load initial condition
        %f = f0;
                
        % Initialize vector variables
        %u_next = zeros(ny,nx);
        
        for tsteps = time
            % for visualization
            rr(:,:) = r(1,1,:,:); 
            h = surf(rr);
        
            % compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_2d(r,u1,u2,v1,v2,t,theta);
            %f = f_eq;
                              
            % (this part can, and should be done in parallel!)
            for  j = 1:nv
                for i = 1:nv
                % Load our physical domain
                u_eq(:,:) = f_eq(j,i,:,:);
                u(:,:) = f(j,i,:,:);
                                    
                % Compute Theta (smoothness coeficients)
                [rx,ry] = theta2d(u,a(j,i),b(j,i));

                % Compute flux Limiter 
                [phix,phiy] = fluxlimiter2d(rx,ry,1); % using limiter = 1

                % Compute TVD Fluxes:
                [F_l,F_r,G_l,G_r] = TVDflux2d(u,a(j,i),b(j,i),dtdx,dtdy,phix,phiy);

                % Compute next time step
                u_next = u - dtdx*(F_r - F_l) - dtdy*(G_r - G_l) ...
                    + (dt/r_time)*(u_eq-u);

                % BC
                u_next(:,1)  = u_next(:,2);
                u_next(:,nx) = u_next(:,nx-1);
                u_next(1,:)  = u_next(2,:);
                u_next(ny,:) = u_next(ny-1,:);

                % UPDATE info
                u = u_next;
                
                % Going back to f
                f(j,i,:,:) = u(:,:);
                end
            end
            
            % Compute macroscopic moments
            [n,j_x,j_y,E] = macromoments2d(k,w1,w2,f,v1,v2);
            
            % UPDATE macroscopic properties 
            % (here lies a parallel computing challenge)
            [r,u1,u2,t,p] = macroproperties2d(n,j_x,j_y,E,nx,ny,nv,theta);
                        
            %fprintf('%2.4f\n',tsteps)
            drawnow
        end
        
    case{2} % WENO k = 3 i.e. O(h^5) 
            % under construction ...
    otherwise
        error('Order must be between 1 and 2');
end

%% Plot results
%plot(x,u);