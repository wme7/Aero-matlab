%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D Semi-classical Boltzmann-ESBGK Equation to recover Quantum Euler sol.
%
%  $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
%  \cdot\nabla_{\vec x} f =\widehat{\Omega}(f)= -\frac{f-f^{eq}}{\tau}$$
%
%              coded by Manuel Diaz, NTU, 2013.02.25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% 1. Yang, J.Y.; Numerical solutions of ideal quantum gas dynamical flows
%    governed by semiclassical ellipsoidal-statistical distribution. DOI:
%    10.1098/rspa.2013.0413 
% 2. H.T. Huynh; A flux reconstruction approach to high-order schemes
%    including Discontinuous Galerkin methods. AIAA 2007.
% 3. J.S. Hestaven & T. Warburton; Nodal Discontinuous Galerkin Methods,
%    Algorithms, Analysis and Applications. Springer 2008.
% 4. Cockburn & Shu; TVB Runge-Kutta Local Projection Discontinuous
%    Galerkin Finite Element Method for conservation laws II: General
%    Framework. 1989.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear all;  close all; clc;

%% Parameters
    parameter.name      = 'SBBGK1d';% Simulation Name
    parameter.CFL       = 15/100;   % CFL condition
    parameter.f_case    = 1;        % {1}Relaxation model, {2}Euler Limit
    parameter.r_time    = 1/1e5;    % Relaxation time (if f_case = 1)
    parameter.tEnd  	= 0.05;     % End time - Parameter part of ICs
    parameter.theta 	= 0;        % {-1}BE, {0}MB, {1}FD.
    parameter.feq_model = 1;        % {1}UU.  {2}ES.
    parameter.v_discrt  = 1;        % {1}NC,  {2}GH
    parameter.method 	= 1;        % {1}NDG, {2}MDG, {3}CPR,
    parameter.IC        = 7;        % IC:{1}~{14}. See SBBGK_IC1d.m
    parameter.P         = 3;        % Polinomial Degree
    parameter.K         = 8;      % Number of elements
    parameter.RK_degree	= 3;        % Number of RK stages
    parameter.RK_stages	= 3;        % Number of RK stages

% Ploting conditions
    plot_figs = 0;	% 0: no, 1: yes please!
    write_ans = 0;	% 0: no, 1: yes please!  

% ID filename for results
    [ID, IDn] = ID_name(parameter);

% Open a Files to store the Results
if write_ans == 1
    file = fopen(IDn,'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
    fprintf(file, 'TITLE = "%s"\n',ID);
    fprintf(file, 'VARIABLES = "x" "density" "velocity" "energy" "pressure" "temperature" "fugacity"\n');
end

%% Preprocessing
% Domain discretization for Spectral methods
    xgrid = mesh1d([0 1],parameter.K,'LGL',parameter.P);
    x  = xgrid.nodeCoordinates;
    dx = xgrid.elementSize; J = xgrid.Jacobian;
    xc = xgrid.elementCenter; wx = xgrid.weights;


%% Microscopic Velocity Discretization (using Discrete Ordinate Method)
% that is to make coincide discrete values of microscopic velocities with
% values as the value points for using a quadrature method, so that we can
% integrate the velocity probability distribution to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature:
    Vv  = [-20,20]; % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,wv,k] = cotes_xw(Vv(1),Vv(2),nv,5); % Using Netwon Cotes Degree 5
        
    case{2} % Gauss Hermite Quadrature:
    nv = 80;         % nodes desired (the actual value)
    [v,wv] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;           % quadrature constant.
    wv = wv.*exp(v.^2);% weighting function of the Gauss-Hermite quadrature
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual nv value will be computed using 'length' vector function:
    nv = length(v);
% Using D.O.M.
    v = reshape(v,[1,1,nv]);    wv = reshape(wv,[1,1,nv]);
    v = repmat(v,Pp,K);         wv = repmat(wv,Pp,K);

% IC in physical Space
% Semiclassical ICs: Fugacity[z], Velocity[u] and Temperature[t] 
    [z0,u0,t0,p0,rho0,E0,~,~] = SBBGK_IC1d(x,parameter.IC);
    
% Applying Discrete Ordinate Method on ICs:
    [z,ux,t] = apply_DG_DOM(z0,u0,t0,nv);  % Semi-classical IC
    [p,rho,E] = apply_DG_DOM(p0,rho0,E0,nv); % Classical IC

% IC for PDF 'f' by assuming the equilibrium state of the macroscopic IC.
% Using the semiclassical Equilibrium distribuition function:
switch fmodel 
    case{1} % U.U.
        f0 = f_equilibrium_1d(z,ux,v,t,theta);
    case{2} % E.S.
        f0 = f_ES_equilibrium_1d(z,p,rho,ux,v,t,theta);
otherwise 
        error('Order must be between 1 and 2');
end
    
% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1
   figure(1)
   f_temp(:,:) = f0(1,:,:);
   surf(f_temp'); grid on; %set(gca,'xDir','reverse');
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end
% Compute Initial Macroscopic Momemts:
    %[rho,rhoux,E,Wxx] = macromoments_DG_1d(k,w,f0,v,ux); %Just for testing
    %[~,~,~,p] = macroproperties_DG_1d(rho,rhou,E,nx,nv,theta);
    
%% Solver
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.

% Prepare for adaptive time stepping
mindx = min(x(2,:)-x(1,:));
dt = mindx*CFL/max(abs(v(:,1))); 
dtdx = dt/mindx;  % precomputed to save someflops

% Time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the classic
% transport Boltzmann equation will resemble to a pure advection equation.
% Thus WENO, TVD, DG or CPR can be used easyly to compute evolution of the
% information inside the domain: 
tic
switch method
    
    case{5} % Nodal DG
        % Using discrete ordinate method (discrete and constant velocity
        % values in phase-space domain)
        a(:,1) = v(1,1,:);
        
        % Load initial condition
        f = f0;
                        
        for tsteps = time
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1
            % Plot f distribution
            figure(1)
            f_temp(:,:) = f(1,:,:);
            surf(f_temp'); grid on; set(gca,'xDir','reverse');
            xlabel('x - Spatial Domain');
            ylabel('v - Velocity Space');
            zlabel('f - Probability');
            % Plot Macroscopic variables
            figure(2)
            subplot(2,3,1); plot(x,rho(:,:,1),'.'); axis tight; title('Density')
            subplot(2,3,2); plot(x,ux(:,:,1),'.'); axis tight; title('x-Velocity')
            subplot(2,3,3); plot(x,p(:,:,1),'.'); axis tight; title('Pressure')
            subplot(2,3,4); plot(x,z(:,:,1),'.'); axis tight; title('Fugacity')
            subplot(2,3,5); plot(x,t(:,:,1),'.'); axis tight; title('Temperature')
            subplot(2,3,6); plot(x,E(:,:,1),'.'); axis tight; title('Energy')
            end
            % Write Results
            if write_ans == 1 && (mod(tsteps,5*dt) == 0 || tsteps == time(end))
                fprintf(file, 'ZONE T = "time %0.4f"\n', tsteps);
                fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n', Pp*K);
                for i = 1:K % Elements
                    for j = 1:N+1 % Element's nodes
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(j,i),rho(j,i,1),ux(j,i,1),E(j,i,1),p(j,i,1),t(j,i,1),z(j,i,1));
                    end
                end
            end
            
            % Compute equilibrium distribution for the current t_step
            switch fmodel
                case{1} % U.U.
                    f_eq = f_equilibrium_1d(z,ux,v,t,theta);
                case{2} % E.S.
                    f_eq = f_ES_equilibrium_1d(z,p,rho,ux,v,t,theta);
                otherwise
                    error('Order must be between 1 and 2');
            end
                        
            % initialize variables
            u_next = zeros(Pp,K);
            u_eq = zeros(Pp,K);
            u = zeros(Pp,K);
                              
            % (this part can, and should be done in parallel!)
            for i = 1:nv
                % load subcase
                switch f_case
                    case{1} % Relaxation Scheme
                        u_eq(:,:) = f_eq(:,:,i);
                        u(:,:) = f(:,:,i);
                    case{2} % Euler Limit
                        u_eq(:,:) = f_eq(:,:,i);
                        u(:,:) = f_eq(:,:,i);
                end
                
                % Limit initial solution
                u = SlopeLimitN(u);
                
                % 3rd order SSP Runge-Kutta
                % SSP RK Stage 1.
                [rhsu]  = SBBGK_RHS1D(a(i),u,u_eq,r_time);
                u1  = u  + dt*rhsu;
                
                % Limit fields
                u1  = SlopeLimitN(u1);
                
                % SSP RK Stage 2.
                [rhsu]  = SBBGK_RHS1D(a(i),u1,u_eq,r_time);
                u2   = (3*u  + u1  + dt*rhsu )/4;
                
                % Limit fields
                u2  = SlopeLimitN(u2);
                
                % SSP RK Stage 3.
                [rhsu]  = SBBGK_RHS1D(a(i),u2,u_eq,r_time);
                u  = (u  + 2*u2  + 2*dt*rhsu )/3;
                
                % Limit solution
                u_next =SlopeLimitN(u);

                % UPDATE info
                u = u_next;
                
                % Going back to f
                f(:,:,i) = u(:,:);
            end
            
            % Compute macroscopic moments
            [rho,rhoux,E,Wxx] = macromoments_DG_1d(k,w,f,v,ux);
            
            % UPDATE macroscopic properties 
            try
            % (here lies a paralellizing computing challenge)
            [z,ux,t,p] = macroproperties_DG_1d(rho,rhoux,E,Wxx,K,Pp,theta,fmodel);
            catch ME
                ME %#ok<NOPTS> % show error message
                break
            end
                                    
            % Apply DOM
            [z,ux,t] = apply_DG_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DG_DOM(p,rho,E,nv); % Classical variables
            
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
    subplot(2,3,1); plot(x,rho(:,:,1),'-o'); axis tight; title('Density')
    subplot(2,3,2); plot(x,ux(:,:,1),'-o'); axis tight; title('x-Velocity')
    subplot(2,3,3); plot(x,p(:,:,1),'-o'); axis tight; title('Pressure')
    subplot(2,3,4); plot(x,z(:,:,1),'-o'); axis tight; title('Fugacity')
    subplot(2,3,5); plot(x,t(:,:,1),'-o'); axis tight; title('Temperature')
    subplot(2,3,6); plot(x,E(:,:,1),'-o'); axis tight; title('Energy')
end