%% 1D Boltzmann Equation
% Finite Difference solution of the Semi-classical Boltzmann Equation to
% recover Euler macroscopic continuum solution. By Manuel Diaz 2013.1.16
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
clear all;  close all; %clc;

%% Simulation Parameters
name        ='SBBGK1d'; % Simulation Name
CFL         = 0.01;     % CFL condition
r_time      = 1/10000;  % Relaxation time
tEnd        = 0.1;      % End time
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 3;        % for DOM-NC = 1, DOM-GH = 2, DDOM-3pGH = 3
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
nx  = 100;                      % Desided number of points in our domain
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

    case{1} % Newton Cotes Quadrature for D.O.M.
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
    nv = length(v);     v = repmat(v,1,nx);     w = repmat(w,1,nx);
    
    case{2} % Gauss Hermite Quadrature for D.O.M.
    nv = 80;        % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    w = w.*exp(v.^2); % weighting function of the Gauss-Hermite quadrature
    v = repmat(v,1,nx);     w = repmat(w,1,nx);
    
    case{3} % Gauss Hermite Quadrature for Dynamic-D.O.M.
    nv = 3;         % nodes required = 3 (the actual value)
    [c_star,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    w = w.*exp(c_star.^2); % weighting function of the Gauss-Hermite quadrature
    J = sqrt(t0);   % Jacobian for every point in x domain, 'J(x)'.
    v = c_star*J + ones(3,1)*u0;  % transformation: v = a*(C*) + ux
    c_star = repmat(c_star,1,nx); w = repmat(w,1,nx); J = repmat(J,nv,1);
        
    otherwise
        error('Order must be between 1, 2 and 3');
end
    
%% Applying Discrete Ordinate Method on ICs:
[z,ux,t] = apply_DOM(z0,u0,t0,nv);  % Semi-classical IC
[p,~,~] = apply_DOM(p0,rho0,E0,nv); % Classical IC
   
%% Initial Equilibrium Distribution 'f0'
% Compute distribution IC: 'f0' of our mesoscopic method by assuming the
% equilibrium state of the macroscopic IC. We use then the semiclassical
% Equilibrium distribuition function:
switch quad 
    case{1,2} % DOM-NC or DOM-GH
        f0 = f_equilibrium_1d(z,ux,v,t,theta);
    case{3}   % DDOM with 3 points GH
        f0 = f_equilibrium_star_1d(z,c_star,theta);
otherwise 
        error('Order must be between 1, 2 and 3');
end

% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1 && (quad == 1 || quad == 2)
   figure(1)
   surf(f0); grid on;
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end

% Compute Initial Macroscopic Momemts:
switch quad 
    case{1,2} % DOM-NC or DOM-GH
        [rho,rhoux,E] = macromoments1d(k,w,f0,v);
    case{3}   % DDOM with 3 GH points
        [rho,rhoux,E] = macromoments_star_1d(J,k,w,f0,v);
otherwise 
        error('Order must be between 1, 2 and 3');
end
    
%% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.
dt = dx*CFL/max(v(:,1)); 
dtdx = dt/dx;  % precomputed to save some flops

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
        a = v;
        
        % Load initial condition
        f = f0;
                        
        for tsteps = time
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1 && (quad == 1 || quad == 2)
            % Plot f distribution
            figure(1)
            surf(f); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            end
            if plot_figs == 1 
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
                        x(1,i),n(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),z(1,i));
                end
            end
            % compute equilibrium distribution for the current t_step
            switch quad
                case{1,2} % DOM-NC or DOM-GH
                    f_eq = f_equilibrium_1d(z,ux,v,t,theta);
                case{3}   % DDOM with 3 points GH
                    f_eq = f_equilibrium_star_1d(z,c_star,theta);
                otherwise
                    error('Order must be between 1, 2 and 3');
            end
            
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
                [r] = theta1d(u,a(i,:));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(u,a(i,:),dtdx,phi);

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
            switch quad
                case{1,2} % DOM-NC or DOM-GH
                    [rho,rhoux,E] = macromoments1d(k,w,f,v);
                case{3}   % DDOM with 3 GH points
                    [rho,rhoux,E] = macromoments_star_1d(J,k,w,f,c_star);
                otherwise
                    error('Order must be between 1, 2 and 3');
            end
            
            % UPDATE macroscopic properties 
            % (here lies a paralellizing computing chalenge)
            [z,ux,t,p] = macroproperties1d(rho,rhoux,E,nx,nv,theta);
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,~,~] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % compute new v and J values if DDOM is used
            if quad == 3 % if DDOM
                J = sqrt(t(1,:));   % Jacobian
                v = ones(3,1)*J.*c_star + ux; % transformation: v = a*(C*) + ux
            end
            % update drawing
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