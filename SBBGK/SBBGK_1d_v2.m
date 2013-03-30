%% 1D Semi-classical Boltzmann-BGK Equation
% Numerical solution of the Boltzmann-BGK Equation to recover Euler macroscopic
% continuum solution. Coded by Manuel Diaz 2012.10.06
%
% Semi-classical Boltzmann-BGK Transport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f) = - \frac{f-f^{eq}}{\tau}$$
%
clear all; close all; %clc;

%% Simulation Parameters
    name	='SBBGK1d'; % Simulation Name
    CFL     = 10/100;   % CFL condition <- Part of IC's parameters
    f_case  = 1;        % {1} Relaxation Model, {2} Euler Limit
    r_time  = 1/10000;  % Relaxation time
    %tEnd  	= 0.10;     % End time <- Part of IC's parameters
    theta 	= 0;        % {-1} BE, {0} MB, {1} FD.
    fmodel  = 2;        % {1} UU. model, {2} ES model.
    quad   	= 2;        % {1} 200NC , {2} 80GH
    method 	= 3;        % {1} Upwind, {2} TVD, {3} WENO3, {4} WENO5
    IC_case	= 7;        % IC: {1}~{14}. See Euler_IC1d.m
  plot_figs = 1;        % 0: no, 1: yes please!
  write_ans = 0;        % 0: no, 1: yes please!
% Using DG
    P_deg	= 1;        % Polinomial Degree
    Pp      = P_deg+1;  % Polinomials Points
% Using RK integration time step
  RK_stages	= 1;        % Number of RK stages

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
% Semiclassical ICs: Fugacity[z], Velocity[u] and Temperature[t] 
    [z0,u0,t0,p0,rho0,E0,tEnd,~] = SSBGK_IC1d(x,IC_case);

%% Microscopic Velocity Discretization (For DOM)
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
    nv = length(v); % actual number of discrete velocity points
    v = repmat(v,1,nx);     w = repmat(w,1,nx);   % using DOM
    
    case{2} % Gauss Hermite Quadrature:
    nv = 80;          % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;            % quadrature constant.
    w = w.*exp(v.^2); % weighting function of the Gauss-Hermite quadrature
    nv = length(v); % actual number of discrete velocity points
    v = repmat(v,1,nx);     w = repmat(w,1,nx);   % using DOM
    
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
        f0 = f_ES_equilibrium_1d(z,p,rho,ux,v,t,theta);
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
    %[rho,rhoux,E,ne,Wxx] = macromoments1d(k,w,f0,v,ux); %Just for testing
    %[~,~,~,p] = macroproperties1d(rho,rhoux,E,ne,Wxx,nx,theta,fmodel);
    
%% Marching Scheme
% First we need to define how big is our time step. Due to the discrete
% ordinate method the problem is similar to evolve the same problem for
% every mesoscopic velocity.
dt = dx*CFL/max(abs(v(:,1))); 
dtdx = dt/dx;  % precomputed to save someflops

% Time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the classic
% transport Boltzmann equation will resemble to a pure advection equation.
% Thus WENO, TVD, DG or CPR can be used easyly to compute evolution of the
% information inside the domain: 
tic

% Load IC
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
        subplot(2,3,2); plot(x,ux(1,:),'.'); axis tight; title('x-Velocity')
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
                x(1,i),rho(1,i),ux(1,i),E(1,i),p(1,i),t(1,i),z(1,i));
        end
    end
    
    % Compute equilibrium distribution for the current t_step
    f_eq = f_equilibrium_1d(z,ux,v,t,theta);
    
    % preallocate
    f_next = zeros(size(f_eq));
    
    for j = 1:nv
        % Load subcase,
        u_eq = f_eq(j,:); u = f(j,:);
        
        % Collision step,
        u_star = collision(u,u_eq,r_time,dt);
        
        % Advection step,
        u_next = stream(u_star,v(j,:),dtdx,method);
              
        % going back to f,
        f_next(j,:) = u_next;
    end
    
    % Update Information
    f = f_next;
    
    % Compute macroscopic moments
    [rho,rhoux,E,ne,Wxx] = macromoments1d(k,w,f,v,ux);
    
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
    
    % Update figures
    drawnow
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
    subplot(2,3,2); plot(x,ux(1,:),'o'); axis tight; title('x-Velocity')
    subplot(2,3,3); plot(x,p(1,:),'o'); axis tight; title('Pressure')
    subplot(2,3,4); plot(x,z(1,:),'o'); axis tight; title('Fugacity')
    subplot(2,3,5); plot(x,t(1,:),'o'); axis tight; title('Temperature')
    subplot(2,3,6); plot(x,E(1,:),'o'); axis tight; title('Energy')
end