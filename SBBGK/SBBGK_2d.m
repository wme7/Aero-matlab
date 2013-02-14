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
    name   ='SBBGK2d';  % Simulation Name
    CFL    = 50/100;    % CFL condition
    r_time = 1/10000;   % Relaxation time
    tEnd   = 0.20;      % Out time
    theta  = 0;         % {-1} BE, {0} MB, {1} FD.
   fmodel  = 1;         % {1} UU. model, {2} ES model.
    quad   = 2;         % {1} NC , {2} GH
    method = 1;         % {1} TVD, {2} WENO3, {3} WENO5
   IC_case = 5;         % Reimann cases: 1~17
 plot_figs = 0;         % 0: no, 1: yes please!
 write_ans = 1;         % 0: no, 1: yes please!
 % Using DG
    P_deg  = 0;         % Polinomial Degree
    Pp 	   = P_deg+1;   % Polinomials Points
% Using RK integration time step
RK_stages   = 4;        % Number of RK stages

%% Space Discretization
nx = 80; ny = 80;
[X,dx,Y,dy] = grid2d(0,1,nx,0,1,ny);
[x,y] = meshgrid(X,Y);

%% Define a ID name for results file
[ID, IDn] = ID_name(name,theta,nx,P_deg,RK_stages,r_time,IC_case);

%% Open a Files to store the Results
if write_ans == 1
    file = fopen(IDn,'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
    fprintf(file, 'TITLE = "%s"\n',ID);
    fprintf(file, 'VARIABLES = "x" "y" "n" "E" "p" "t" "z"\n');
end

%% Velocity Discretization:
switch quad

    case{1} % Newton Cotes Quadrature:
    V  = [-10,10];  % range: a to b
    nv = 50;        % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % cotes Degree 5

    case{2} % Gauss Hermite Quadrature:
    nv = 20;        % nodes desired (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
     k = 1;         % quadrature constant.
     w = w.*exp(v.^2);
    
    otherwise
        error('Order must be between 1 and 2');
end
%% Velocity-Space Grid:
% The actual nv value will be computed using 'length' vector function:
nv = length(v); 

% Using Discrete Ordinate Method (Replicating Data)
    [vx,vy] = meshgrid(v,v');       [wx,wy] = meshgrid(w,w');
    vx = repmat(vx,[1,1,ny,nx]);    wx = repmat(wx,[1,1,ny,nx]);
    vy = repmat(vy,[1,1,ny,nx]);    wy = repmat(wy,[1,1,ny,nx]);

%% Initial Conditions
% Load Macroscopic Velocity, Temperature and Fugacity
    [z0,ux0,uy0,t0] = SSBGK_IC2d(x,y,IC_case);
    
% Using Discrete Ordinate Method: (Replicating Data for easyer arrays Ops)
    z = zeros(nv,nv);   ux = zeros(nv,nv);
    uy = zeros(nv,nv);   t = zeros(nv,nv);

for j = 1:ny
     for i = 1:nx
     z(:,:,j,i) = z0(j,i)*ones(nv,nv);  ux(:,:,j,i) = ux0(j,i)*ones(nv,nv);
     uy(:,:,j,i) = uy0(j,i)*ones(nv,nv);  t(:,:,j,i) = t0(j,i)*ones(nv,nv);
     end
end

% Compute distribution IC of our mesoscopic method by assuming the equilibrium
% state of the macroscopic IC. Using the semiclassical Equilibrium
% distribuition function:
    f0 = f_equilibrium_2d(z,ux,uy,t,vx,vy,theta);

% Plot Probability Distribution Function, f, for an expecific point in space:
if plot_figs == 1
    figure(1)  
    i = 2; j = 2;
    zz = f0(:,:,j,i);
    surf(zz); grid on;
    xlabel('v1 - Velocity Space in x'); 
    ylabel('v2 - Velocity Space in y');
	zlabel('f - Probability');
	clear zz i j
end
    
% Compute Initial Macroscopic Momemts:
    [rho,rhoux,rhouy,E] = macromoments2d(k,wx,wy,f0,vx,vy);
    [~,~,~,~,p] = macroproperties2d(rho,rhoux,rhouy,E,nx,ny,theta);
    
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
tic
switch method
    
    case{1} % TVD 0(h^2)
        % Using Discrete Ordinate Method (Discrete and constant velocity
        % values in phase-space domain)
        a = vx(:,:,1,1);
        b = vy(:,:,1,1);
        
        % Load initial condition
        f = f0;
                
        % Initialize vector variables
        u_next = zeros(ny,nx);
        u_eq = zeros(ny,nx);
        u = zeros(ny,nx);
        
        for tsteps = time
            % Plot and redraw figures every time step for visualization
            if plot_figs == 1
            % Plot Macroscopic variables
            figure(2)
            subplot(2,3,1); nn(:,:) = rho(:,:); h2 = surf(nn); title('Density')
            subplot(2,3,2); pp(:,:) = p(:,:); h4 = surf(pp); title('Pressure')
            subplot(2,3,3); tt(:,:) = t(1,1,:,:); h3 = surf(tt); title('Temperature')
            subplot(2,3,4); zz(:,:) = z(1,1,:,:); h1 = surf(zz); title('Fugacity')
            subplot(2,3,5); EE(:,:) = E(:,:); h5 = surf(EE); title('Energy')
            end
            % Write Results
            if write_ans == 1 && (mod(tsteps,10*dt) == 0 || tsteps == time(end))
                fprintf(file, 'ZONE T = "time %0.4f"\n', tsteps);
                fprintf(file, 'I = %d, J = %d, K = 1, F = POINT\n\n',nx,ny);
                for j = 1:ny
                for i = 1:nx
                    fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
                        x(j,i),y(j,i),rho(j,i),E(j,i),p(j,i),t(1,1,j,i),z(1,1,j,i));
                end
                end
            end
            
            % Compute Equilibrium Distribution for the current time step
            f_eq = f_equilibrium_2d(z,ux,uy,t,vx,vy,theta);
                                          
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
                    - (dt/r_time)*(u-u_eq);

                % BCs
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
            [rho,rhoux,rhouy,E] = macromoments2d(k,wx,wy,f,vx,vy);
            
            % UPDATE macroscopic properties 
            % (here lies a parallel computing challenge)
            [z1,ux1,uy1,t1,p1] = macroproperties2d(rho,rhoux,rhouy,E,nx,ny,theta);
            
            % Using Discrete Ordinate Method: (Replicating data for easier arrays Ops)
            for j = 1:ny
                for i = 1:nx
                z(:,:,j,i)  = z1(j,i)*ones(nv,nv);  ux(:,:,j,i) = ux1(j,i)*ones(nv,nv);
                uy(:,:,j,i) = uy1(j,i)*ones(nv,nv); t(:,:,j,i) = t1(j,i)*ones(nv,nv);
                end
            end
            
            % Update Figures
            drawnow
        end
        
    case{2} % WENO(r=2)3, O(h^5) 
            % under construction ...
    otherwise
        error('Order must be between 1 and 2');
end
toc

%% Write Results
fprintf('Simulation has been completed succesfully!\n')
if write_ans == 1
    fclose(file);
    fprintf('All Results have been saved!\n')
end

if plot_figs ~= 1
    % Plot Macroscopic variables
    figure(2)
    subplot(2,3,1); nn(:,:) = rho(:,:); h2 = surface(nn); title('Density')
    subplot(2,3,2); pp(:,:) = p(:,:); h4 = surface(pp); title('Pressure')
    subplot(2,3,3); tt(:,:) = t(1,1,:,:); h3 = surface(tt); title('Temperature')
    subplot(2,3,4); zz(:,:) = z(1,1,:,:); h1 = surface(zz); title('Fugacity')
    subplot(2,3,5); EE(:,:) = E(:,:); h5 = surface(EE); title('Energy')
end