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
CFL         = 5/100;    % CFL condition <- can be made part of IC's parameters
f_case      = 1;        % {1} Relaxation Model, {2} Euler limit
r_time      = 1/10000;  % Relaxation time
%tEnd       = 0.04;     % End time <- part of IC's parameters
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 3;        % {1} DOM-200NC, {2} DOM-80GH, {3} DDOM-3GH
method      = 1;        % {1} Upwind, {2} TVD, {3} WENO3, {4} WENO5
fmodel      = 1;        % {1} UU model, {2} ES model
IC_case     = 1;        % % IC: {1}~{14}. See Euler_IC1d.m
plot_figs   = 1;        % 0: no, 1: yes please!
write_ans   = 0;        % 0: no, 1: yes please!
% If using DG
P_deg       = 1;        % Polinomial Degree
Pp          = P_deg+1;  % Polinomials Points
% Using RK integration time step
RK_stages   = 1;        % Number of RK stages

%% Space Discretization
nx  = 100;                      % Desided number of points in our domain
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Define a ID name for results file
[ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,r_time,IC_case,fmodel,f_case,method);

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
    [z0,u0,t0,p0,rho0,E0,tEnd,~] = SSBGK_IC1d(x,IC_case);
   
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
    J = ones(size(v)); % Jacobian = 1;
    
    case{2} % Gauss Hermite Quadrature for D.O.M.
    nv = 80;        % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    w = w.*exp(v.^2); % weighting function of the Gauss-Hermite quadrature
    v = repmat(v,1,nx);     w = repmat(w,1,nx);
    J = ones(size(v)); % Jacobian = 1;
    
    case{3} % Gauss Hermite Quadrature for Dynamic-D.O.M.
    nv = 3;         % nodes required = 3 (the actual value)
    [c_star,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    w = w.*exp(c_star.^2); % weighting function of the Gauss-Hermite quadrature
    J = sqrt(t0);   % Jacobian for every point in x domain, 'J(x)'.
    v = c_star*J + ones(nv,1)*u0;  % transformation: v = a*(C*) + ux
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
    f0 = f_equilibrium_1d(z,ux,v,t,theta);

        % Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1 %&& (quad == 1 || quad == 2)
   figure(1)
   surf(f0); grid on;
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end

% Compute Initial Macroscopic Momemts:
[rho,rhoux,E,ne] = macromoments_star_1d(J,k,w,f0,v,ux);
    
%% Marching Scheme

% Initial time
time = 0;

% Initial time step 'dt'
dt = dx*CFL/max(abs(v(:,1)));

% Set 50 save points
savept = linspace(0,tEnd,50);
 
tic
switch method
            
    case{1} % Upwind O(h)
        % Load IC
        f = f0;
        
        % Main loop
        %for tsteps = time
        while time <= tEnd
            % Update discrete velocity points for DDOM
            a = v;
            
            % Update time step
            dt = dx*CFL/max(max(abs(a)));
            time = time + dt;
            dtdx = dt/dx; 
            
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1 
            % Plot f distribution
            figure(1)
            surf(f); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            end
            if plot_figs == 1 
            % Plot Macroscopic variables
            figure(2) %range = ([0,1,0,1.4])
            subplot(2,3,1); plot(x,rho(1,:),'.'); axis auto; title('Density')
            subplot(2,3,2); plot(x,ux(1,:),'.'); axis auto; title('velocity in x')
            subplot(2,3,3); plot(x,p(1,:),'.'); axis auto; title('Pressure')
            subplot(2,3,4); plot(x,z(1,:),'.'); axis auto; title('Fugacity')
            subplot(2,3,5); plot(x,t(1,:),'.'); axis auto; title('Temperature')
            subplot(2,3,6); plot(x,E(1,:),'.'); axis auto; title('Energy')
            end
            % Write Results
            tol = dt/2;
            %if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
            if write_ans == 1 && max( (time < savept+tol) & (time > savept-tol) )
                fprintf(file, 'ZONE T = "time %0.4f"\n', time);
                fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);
                for i = 1:nx
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(1,i),rho(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),z(1,i));
                end
            end
            % Compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(z,ux,v,t,theta);
            
            % initialize variables
            u_eq = zeros(1,nx);
            u  = zeros(1,nx);
                            
            % (this part can be done in parallel!)
            for i = 1:nv       
                % Load subcase
                switch f_case
                    case{1} % Relaxation Scheme
                        u_eq(:) = f_eq(i,:);
                        u(:) = f(i,:);
                    case{2} % Euler Limit
                        u_eq(:) = f_eq(i,:);
                        u(:) = f_eq(i,:);
                end
                
                % Compute Upwind Fluxes
                [F_left,F_right] = Upwindflux1d(u,a(i,:));
             
                % Compute next time step
                 u_next = u - dtdx*(F_right - F_left) ...
                     - (dt/r_time)*(u-u_eq);
                % BC
                u_next(1) = u_next(2);
                u_next(nx) = u_next(nx-1);

                % UPDATE info
                u = u_next;
                             
                % Going back to f
                f(i,:) = u(:);
            end
            
            % Compute macroscopic moments
            [rho,rhoux,E,ne,Wxx] = macromoments_star_1d(J,k,w,f,v,ux);
            
            % UPDATE macroscopic properties 
            try
            % (here lies a paralellizing computing challenge)
            [z,ux,t,p] = macroproperties1d(rho,rhoux,E,ne,Wxx,nx,theta,fmodel);
            catch ME
                ME %#ok<NOPTS> % show error message
                break
            end
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Compute new v and J values if DDOM is used
            if quad == 3 % if DDOM
                J = sqrt(t);   % Jacobian
                v = J.*c_star + ux; % transformation: v = a*(C*) + ux
            end
            % Update figures
            drawnow
        end
        
    case{2} % TVD ... coming soon
        % Load IC
        f = f0;
        
        % Main loop
        for tsteps = time
            % Update discrete velocity points for DDOM
            a = v;
            
            % Update time step
            dt = dx*CFL/max(max(abs(a)));
            time = time + dt;
            dtdx = dt/dx;
            
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1 
            % Plot f distribution
            figure(1)
            surf(f); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            end
            if plot_figs == 1 
            % Plot Macroscopic variables
            figure(2) %range = ([0,1,0,1.4])
            subplot(2,3,1); plot(x,rho(1,:),'.'); axis auto; title('Density')
            subplot(2,3,2); plot(x,ux(1,:),'.'); axis auto; title('velocity in x')
            subplot(2,3,3); plot(x,p(1,:),'.'); axis auto; title('Pressure')
            subplot(2,3,4); plot(x,z(1,:),'.'); axis auto; title('Fugacity')
            subplot(2,3,5); plot(x,t(1,:),'.'); axis auto; title('Temperature')
            subplot(2,3,6); plot(x,E(1,:),'.'); axis auto; title('Energy')
            end
            % Write Results
            %if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
            if write_ans == 1 && max( (time < savept+tol) & (time > savept-tol) )
                fprintf(file, 'ZONE T = "time %0.4f"\n', time);
                fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);
                for i = 1:nx
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(1,i),rho(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),z(1,i));
                end
            end
            % Compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(z,ux,v,t,theta);
            
            % initialize variables
            u_eq = zeros(1,nx);
            u  = zeros(1,nx);
                            
            % (this part can be done in parallel!)
            for i = 1:nv       
                % Load subcase
                switch f_case
                    case{1} % Relaxation Scheme
                        u_eq(:) = f_eq(i,:);
                        u(:) = f(i,:);
                    case{2} % Euler Limit
                        u_eq(:) = f_eq(i,:);
                        u(:) = f_eq(i,:);
                end
                
                % Compute the smoothness factors, r(j), from data, u(j).
                [r] = theta1d(u,a(i,:));

                % Compute the Flux Limiter
                [phi] = fluxlimiter1d(r,1); % using limiter = 1

                % Compute TVD Fluxes
                [F_left,F_right] = TVDflux1d(u,a(i,:),dtdx,phi);
             
                % Compute next time step
                 u_next = u - dtdx*(F_right - F_left) ...
                     - (dt/r_time)*(u-u_eq);
                % BC
                u_next(1) = u_next(2);
                u_next(nx) = u_next(nx-1);

                % UPDATE info
                u = u_next;
                             
                % Going back to f
                f(i,:) = u(:);
            end
            
            % Compute macroscopic moments
            [rho,rhoux,E,ne,Wxx] = macromoments_star_1d(J,k,w,f,v,ux);
            
            % UPDATE macroscopic properties 
            try
            % (here lies a paralellizing computing challenge)
            [z,ux,t,p] = macroproperties1d(rho,rhoux,E,ne,Wxx,nx,theta,fmodel);
            catch ME
                ME %#ok<NOPTS> % show error message
                break
            end
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Compute new v and J values if DDOM is used
            if quad == 3 % if DDOM
                J = sqrt(t);   % Jacobian
                v = J.*c_star + ux; % transformation: v = a*(C*) + ux
            end
            % Update figures
            drawnow
        end
        
    case{3} % WENO(r = 2)3, 0(h^5)
        % Load IC
        f = f0;
        
        % Main loop
        for tsteps = time
            % Update discrete velocity points for DDOM
            a = v;
            
            % Update time step
            dt = dx*CFL/max(max(abs(a)));
            time = time + dt;
            dtdx = dt/dx;
            
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1 
            % Plot f distribution
            figure(1)
            surf(f); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            end
            if plot_figs == 1 
            % Plot Macroscopic variables
            figure(2) %range = ([0,1,0,1.4])
            subplot(2,3,1); plot(x,rho(1,:),'.'); axis auto; title('Density')
            subplot(2,3,2); plot(x,ux(1,:),'.'); axis auto; title('velocity in x')
            subplot(2,3,3); plot(x,p(1,:),'.'); axis auto; title('Pressure')
            subplot(2,3,4); plot(x,z(1,:),'.'); axis auto; title('Fugacity')
            subplot(2,3,5); plot(x,t(1,:),'.'); axis auto; title('Temperature')
            subplot(2,3,6); plot(x,E(1,:),'.'); axis auto; title('Energy')
            end
            % Write Results
            %if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
            if write_ans == 1 && max( (time < savept+tol) & (time > savept-tol) )
                fprintf(file, 'ZONE T = "time %0.4f"\n', time);
                fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', nx);
                for i = 1:nx
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(1,i),rho(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),z(1,i));
                end
            end
            % Compute equilibrium distribution for the current t_step
            f_eq = f_equilibrium_1d(z,ux,v,t,theta);
            
            % initialize variables
            u_eq = zeros(1,nx);
            u  = zeros(1,nx);
                            
            % (this part can be done in parallel!)
            for i = 1:nv       
                % Load subcase
                switch f_case
                    case{1} % Relaxation Scheme
                        u_eq(:) = f_eq(i,:);
                        u(:) = f(i,:);
                    case{2} % Euler Limit
                        u_eq(:) = f_eq(i,:);
                        u(:) = f_eq(i,:);
                end
                
                % Scalar Flux Spliting
                [vp,vn] = WENO_scalarfluxsplit(a(i)*u);
                
                % Reconstruct Fluxes values at cells interfaces
                [F_left,F_right] = WENO3_1d_driver(vp,vn);
             
                % Compute next time step
                 u_next = u - dtdx*(F_right - F_left) ...
                     - (dt/r_time)*(u-u_eq);
                % BC
                u_next = WENO3_1d_BCs(u_next,2,nx); %2:Neumann BC

                % UPDATE info
                u = u_next;
                             
                % Going back to f
                f(i,:) = u(:);
            end
            
            % Compute macroscopic moments
            [rho,rhoux,E,ne,Wxx] = macromoments_star_1d(J,k,w,f,v,ux);
            
            % UPDATE macroscopic properties 
            try
            % (here lies a paralellizing computing challenge)
            [z,ux,t,p] = macroproperties1d(rho,rhoux,E,ne,Wxx,nx,theta,fmodel);
            catch ME
                ME %#ok<NOPTS> % show error message
                break
            end
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Compute new v and J values if DDOM is used
            if quad == 3 % if DDOM
                J = sqrt(t);   % Jacobian
                v = J.*c_star + ux; % transformation: v = a*(C*) + ux
            end
            % Update figures
            drawnow
        end
        
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

if plot_figs ~= 1
    % Plot Macroscopic variables
    figure(2)
    subplot(2,3,1); plot(x,rho(1,:),'o'); axis tight; title('Density')
    subplot(2,3,2); plot(x,ux(1,:),'o'); axis tight; title('velocity in x')
    subplot(2,3,3); plot(x,p(1,:),'o'); axis tight; title('Pressure')
    subplot(2,3,4); plot(x,z(1,:),'o'); axis tight; title('Fugacity')
    subplot(2,3,5); plot(x,t(1,:),'o'); axis tight; title('Temperature')
    subplot(2,3,6); plot(x,E(1,:),'o'); axis tight; title('Energy')
end