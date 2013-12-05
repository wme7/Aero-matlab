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
    parameter.v_discrt  = 2;        % {1}NC,  {2}GH
    parameter.method 	= 1;        % {1}NDG, {2}MDG, {3}CPR,
    parameter.IC        = 7;        % IC:{1}~{14}. See SBBGK_IC1d.m
    parameter.P         = 3;        % Polinomial Degree
    parameter.K         = 10;       % Number of elements
    parameter.RK_degree	= 3;        % Number of RK stages
    parameter.RK_stages	= 3;        % Number of RK stages

% Ploting conditions
    plot_figs = 1;	% 0: no, 1: yes please!
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
% Physical domain discretization for Spectral methods
    xgrid = mesh1d([0 1],parameter.K,'LGL',parameter.P);
    x  = xgrid.nodeCoordinates; nx = xgrid.nNodes;
    dx = xgrid.elementSize; J = xgrid.Jacobian;
    xc = xgrid.elementCenter; wx = xgrid.weights;

% DG tools
switch parameter.method
    case 1 % Nodal DG
        tool = DGtools(xgrid.solutionPoints);
        V = tool.Vandermonde2;
        invM = tool.nodalInvMassMatrix;
        Dr = tool.nodalCoefDiffMatrix;
    case 2 % Modal DG
        tool = DGtools(xgrid.solutionPoints);
        V = tool.Vandermonde;
        Emat = zeros(K+1,2);
        Emat(1,1)=1; Emat(K+1,2)=1;
        Lift = V.(V'*Emat);
    case 3 % CPR
        RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
        dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);
        l = LagrangePolynomial(xgrid.solutionPoints);
        L.lcoef = double(subs(l.lagrangePolynomial,-1));
        L.rcoef = double(subs(l.lagrangePolynomial,1));
        L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));
end

% Phase domain discretization for DOM
    BGK = BGKtools(parameter);
    v=BGK.c;    nv=BGK.nc;  k=BGK.kc;   wv=BGK.wc;

% IC for physical space (Macroscopic variables)
    % Classical ICs: Pressure[p], Density[rho] and Energy[E],
    % Semiclassical ICs: Fugacity[z], Velocity[u] and Temperature[t].
    [z0,u0,t0,p0,rho0,E0,~,~] = SBBGK_IC1d(x,parameter.IC);
    
% Applying Discrete Ordinate Method on ICs:
    [z,ux,t] = BGK.apply_DG_DOM(z0,u0,t0,nv);
    [p,rho,E] = BGK.apply_DG_DOM(p0,rho0,E0,nv);

% IC for the PDF 'f':
	% Assume the equilibrium state of the physical IC.
	f0 = BGK.f_equilibrium_1d(z,ux,v,t,parameter.theta);
    
% Plot IC PDF 'f', in Phase-Space:
if plot_figs==1
   figure(1)
   f_temp(:,:) = f0(3,:,:);
   surf(f_temp'); grid on; set(gca,'xDir','reverse');
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end

% Compute initial Macroscopic momemts
    % Just for testing purposes
    [rho,rhoux,E,Wxx] = macromoments_DG_1d(k,wv,f0,v,ux); 
    [~,~,~,p] = macroproperties_DG_1d(rho,rhoux,E,nx,nv,parameter.theta);
    
%% Solver Loop
% First we need to define how big is our time step. Due to the discrete
% ordinate method, the problem is similar to evolve the a similar problem
% for every mesoscopic velocity, 'v'.

% Prepare for adaptive time stepping
dt = dx*CFL/max(abs(v(:,1))); 
dtdx = dt/mindx;  % precomputed to save someflops

% Time domain discretization
time = 0:dt:tEnd;

% By negleting any force field acting over our domian, the transport
% Boltzmann equation will resembles to a pure advection equation. 
% Thus NDG, MDG or CPR can be used easily to compute evolution of the
% information inside the domain.

% DOM vector of discrete velocities
a(:,1) = v(1,1,:);

% Load initial condition
f = f0;

tic
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
    f_eq = f_ES_equilibrium_1d(z,p,rho,ux,v,t,theta);

    % initialize variables
    u_next = zeros(Pp,K);
    u_eq = zeros(Pp,K);
    u = zeros(Pp,K);

    % (this part can, and should be done in parallel!)
    for i = 1:nv
        % load subcase
        switch f_case
            case 1 % Relaxation Scheme
                u_eq(:,:) = f_eq(:,:,i);
                u(:,:) = f(:,:,i);
            case 2 % Euler Limit
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
toc

%% Save file with Results
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