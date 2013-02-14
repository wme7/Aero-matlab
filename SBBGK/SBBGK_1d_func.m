function SBBGK_1d_func(name,CFL,r_time,tEnd,theta,quad,method,IC_case, ...
        plot_figs,write_ans,P_deg,Pp,RK_stages,nx)
%% 1D Boltzmann Equation
% Finite Difference solution of the Semi-classical Boltzmann equation to
% recover Euler macroscopic continuum solution. By Manuel Diaz 2012.10.06
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
%clc;  clear all;  close all;

%% Simulation Parameters
% name        ='SBBGK1d'; % Simulation Name
% CFL         = 15/100;     % CFL condition
% r_time      = 1/10000;  % Relaxation time
% tEnd        = 0.02;      % End time
% theta       = 0;        % {-1} BE, {0} MB, {1} FD.
% quad        = 1;        % for NC = 1 , GH = 2
% method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
% IC_case     = 5;        % Reimann IC cases available :{1,2,3,4,5,6,7}
% plot_figs   = 0;        % 0: no, 1: yes please!
% write_ans   = 0;        % 0: no, 1: yes please!
% % Using DG
% P_deg       = 0;        % Polinomial Degree
% Pp          = P_deg+1;  % Polinomials Points
% % Using RK integration time step
% RK_stages   = 4;        % Number of RK stages

%% Space Discretization
% nx  = 100;                      % Desided number of points in our domain
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Define a ID name for results file
[ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,r_time,IC_case);

%% Open a Files to store the Results
file = fopen(IDn,'w');
% 'file' gets the handel for the file "case.plt".
% 'w' specifies that it will be written.
% similarly 'r' is for reading and 'a' for appending.

fprintf(file, 'TITLE = "%s"\n',ID);
fprintf(file, 'VARIABLES = "x" "density" "velocity" "energy" "pressure" "temperature" "fugacity"\n');

%% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-40,40];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    nv = 80;          % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;            % quadrature constant.
    w = w.*exp(v.^2); % weighting fucntion of the Gauss-Hermite quadrature
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual nv value will be computed using 'lenght' vector function:
nv = length(v); 

% Using D.O.M.
    v = repmat(v,1,nx);     w = repmat(w,1,nx);

% Initialize Arrays
%    ux = zeros(nv,nx);      t = zeros(nv,nx);       
%    r = zeros(nv,nx);       n = zeros(nv,nx);
	p = zeros(nv,nx);   


%% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
    [r0,u0,t0] = reimann_IC1d(x,IC_case);
    
% Using Discrete Ordinate Method:
    r = repmat(r0,nv,1); ux = repmat(u0,nv,1); t = repmat(t0,nv,1);

% Compute distribution IC of our mesoscopic method by assuming the equilibrium 
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
    f0 = f_equilibrium_1d(r,ux,v,t,theta);

% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1
   figure(1)
   surf(f0); grid on;
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end
% Compute Initial Macroscopic Momemts:
    [n,j_x,E] = macromoments1d(k,w,f0,v);
    
%% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.
dt = dx*CFL/max(v(:,1)); 
dtdx = dt/dx;  % precomputed to save someflops

% Time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the classic
% transport Boltzmann equation will resemble to a pure advection equation.
% Thus WENO, TVD, DG or CPR can be used easyly to compute evolution of the
% information inside the domain: 
tic
switch method
    
    case{1} % TVD 0(h^2)
        % Using discrete ordinate method (discrete and constant velocity
        % values in phase-space domain)
        a = v(:,1);
        
        % Load initial condition
        f = f0;
                        
        for tsteps = time
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1
            figure(2)
            subplot(2,3,1); plot(x,n(1,:),'.'); axis tight; title('Density')
            subplot(2,3,2); plot(x,p(1,:),'.'); axis tight; title('Pressure')
            subplot(2,3,3); plot(x,t(1,:),'.'); axis tight; title('Temperature')
            subplot(2,3,4); plot(x,r(1,:),'.'); axis tight; title('Fugacity')
            subplot(2,3,5); plot(x,ux(1,:),'.'); axis tight; title('velocity in x')
            subplot(2,3,6); plot(x,E(1,:),'.'); axis tight; title('Energy')
            end
            % Write Results
            if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
                fprintf(file, 'ZONE T = "time %0.4f"\n', tsteps);
                fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);
                for i = 1:nx
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(1,i),n(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),r(1,i));
                end
            end
            % compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(r,ux,v,t,theta);
            
            % initialize variables
            u_next = zeros(1,nx);
            u_eq = zeros(1,nx);
            u = zeros(1,nx);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                u_eq(:) = f_eq(i,:);
                u(:) = f(i,:);
                
                % Compute the smoothness factors, r(j), from data, u(j).
                [r] = theta1d(u,a(i));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(u,a(i),dtdx,phi);

                % Compute next time step
                u_next = u - dtdx*(F_right - F_left) ...
                    + (dt/r_time)*(u_eq-u);

                % BC
                u_next(1) = u_next(2);
                u_next(nx) = u_next(nx-1);

                % UPDATE info
                u = u_next;
                
                % Going back to f
                f(i,:) = u(:);
            end
            
            % Compute macroscopic moments
            [n,j_x,E] = macromoments1d(k,w,f,v);
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing chalenge)
            [r,ux,t,p] = macroproperties1d(n,j_x,E,nx,nv,theta);
            
            drawnow
        end
        
    case{2} % WENO k = 3 i.e. O(h^5) 
        
    otherwise
        error('Order must be between 1 and 2');
end
toc
%% Close file with Results
fclose(file);
fprintf('Simulation has been completed succesfully!\n')
fprintf('All Results have been saved!\n')