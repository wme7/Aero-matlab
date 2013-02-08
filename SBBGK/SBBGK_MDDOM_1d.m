%% 1D Boltzmann Equation
% Finite Difference solution of the Semi-classical Boltzmann Equation to
% recover Euler macroscopic continuum solution. By Manuel Diaz 2012.10.06
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
clc;  clear all;  close all;

%% Simulation Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.05;     % CFL condition
r_time      = 1/1000;  % Relaxation time
tEnd        = 0.1;      % End time
theta       = 1;        % {-1} BE, {0} MB, {1} FD.
fmodel      = 1;        % {1} UU. model, {2} ES model.
quad        = 2;        % for DOM-NC = 1, DOM-GH = 2, MD-DOM-GH = 3
method      = 1;        % for TVD = 1, WENO3 = 2, WENO5 = 3
IC_case     = 1;        % IC: {1}Sod's, {2}LE, {3}RE, {4}DS, {5}SS, {6}Cavitation
plot_figs   = 1;        % 0: no, 1: yes please!
write_ans   = 0;        % 0: no, 1: yes please!
% Using DG
P_deg       = 0;        % Polinomial Degree
Pp          = P_deg+1;  % Polinomials Points
% Using RK integration time step
RK_stages   = 4;        % Number of RK stages

%% Space Discretization
nx  = 160;                      % Desided number of points in our domain
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Define a ID name for results file
[ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,r_time,IC_case);

%% Open a Files to store the Results
if write_ans == 1
    file = fopen(IDn,'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
    fprintf(file, 'TITLE = "%s"\n',ID);
    fprintf(file, 'VARIABLES = "x" "density" "velocity" "energy" "pressure" "temperature" "fugacity"\n');
end

%% Initial Conditions in physical Space
% Load Macroscopic Fugacity [z], Velocity[u] and Temperature[t] 
    [z0,u0,t0,p0,rho0,E0] = SSBGK_IC1d(x,IC_case);

%% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
    nv = length(v);     v = repmat(v,1,nx);     w = repmat(w,1,nx);
        
    case{2} % Gauss Hermite Quadrature:
    nv = 60;          % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;            % quadrature constant.
    w = w.*exp(v.^2); % weighting function of the Gauss-Hermite quadrature
    v = repmat(v,1,nx);     w = repmat(w,1,nx);
    
    case{3} % Uniformly distributed velocity points
    V  = [-7,7]; % range: a to b
    nv = 20;       % nodes desired (the actual value)
    v = linspace(V(1),V(2),nv)';
    [v_star,w] = GaussHermite(3); % for integrating range: -inf to inf
    w = w.*exp(v_star.^2); % weighting function of the Gauss-Hermite quadrature
    k = 1;         % quadrature constant.
    v = repmat(v,1,nx); v_star = repmat(v_star,1,nx); w = repmat(w,1,nx);
    
    otherwise
        error('Order must be between 1 and 2');
end

%% Applying Discrete Ordinate Method on ICs:
[z,ux,t] = apply_DOM(z0,u0,t0,nv);  % Semi-classical IC
[p,rho,E] = apply_DOM(p0,rho0,E0,nv); % Classical IC

% Compute distribution IC of our mesoscopic method by assuming the equilibrium 
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
switch fmodel 
    case{1} % U.U.
        f0 = f_equilibrium_1d(z,ux,v,t,theta);
    case{2} % E.S.
        f0 = f_SE_equilibrium_1d(z,p,rho,ux,v,t,theta);
otherwise 
        error('Order must be between 1 and 2');
end
    
% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1
   figure(1)
   surf(f0); grid on; set(gca,'xDir','reverse');
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end
% Compute Initial Macroscopic Momemts:
%    [rho,rhou,E] = macromoments1d(k,w,f0,v); Just for testing
    
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
            % Plot f distribution
            figure(1)
            surf(f); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            % Plot Macroscopic variables
            figure(2)
            subplot(2,3,1); plot(x,rho(1,:),'.'); axis tight; title('Density')
            subplot(2,3,2); plot(x,ux(1,:),'.'); axis tight; title('velocity in x')
            subplot(2,3,3); plot(x,p(1,:),'.'); axis tight; title('Pressure')
            subplot(2,3,4); plot(x,z(1,:),'.'); axis tight; title('Fugacity')
            subplot(2,3,5); plot(x,t(1,:),'.'); axis tight; title('Temperature')
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
            
            % Compute equilibrium distribution for the current t_step
            switch fmodel
                case{1} % U.U.
                    f_eq = f_equilibrium_1d(z,ux,v,t,theta);
                case{2} % E.S.
                    f_eq = f_SE_equilibrium_1d(z,p,rho,ux,v,t,theta);
                otherwise
                    error('Order must be between 1 and 2');
            end
                        
            % initialize variables
            u_next = zeros(1,nx);
            u_eq = zeros(1,nx);
            u = zeros(1,nx);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                u_eq(:) = f_eq(i,:);
                u(:) = f_eq(i,:);
                
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
            
            % Compute Macroscopic Moments
            switch quad
                case{1,2} % DOM
                    [rho,rhoux,E] = macromoments1d(k,w,f,v);
                case{3}   % iDOM
                    [rho,rhoux,E] = macromoments1d_iDOM(k,w,f,v,v_star);
                otherwise
                    error('case not available')
            end
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing challenge)
            [z,ux,t,p] = macroproperties1d(rho,rhoux,E,nx,nv,theta);
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Update figures
            drawnow
        end
        
    case{2} % WENO k = 3 i.e. O(h^5) 
        
    otherwise
        error('Order must be between 1 and 2');
end
toc
%% Close file with Results
fprintf('Simulation has been completed succesfully!\n')
if write_ans == 1
    fclose(file);
    fprintf('All Results have been saved!\n')
end