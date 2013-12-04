%% 1D Boltzmann Equation
% Finite Difference solution of the Semi-classical Boltzmann Equation to
% recover Euler macroscopic continuum solution. By Manuel Diaz 2013.3.13
%
% Boltzmann Trasnport Equation:
%
% $$\frac{\partial f}{\partial t}+\vec F\cdot \nabla_p f + \vec v
% \cdot\nabla_{\vec x} f =\widehat{\Omega } (f)$$
%
% Based on ideas of,
% CT Hsu, S.W. Chiang and K.F. Sin
% A Novel Dynamic Quadrature Scheme for Solving Boltzmann Equation with
% DOM. Doi: 10.4208/cicp.150510.150511s
%
clear all;  close all; %clc;

%% Simulation Parameters
CFL         = 10/100;   % CFL condition <- can be made part of IC's
f_case      = 1;        % {1}Relaxation Model, {2}Euler limit
r_time      = 1/10000;  % Relaxation time
%tEnd       = 0.20;     % End time <- part of IC's
theta       = 0;        % {-1} BE, {0} MB, {1} FD.
quad        = 3;        % {1} DOM-200NC, {2} DOM-80GH, {3} DDOM-3GH
method      = 1;        % for {1} Upwind
fmodel      = 1;        % UU model for f^Eq <- fixed for the moment
IC_case     = 7;        % % IC: {1}~{15}. See Euler_IC1d.m
plot_figs   = 1;        % 0: no, 1: yes please!

%% Space Discretization
nx  = 100;                      % Desided number of points in our domain
x   = linspace(0,1,nx);         % Physical domain -x
dx  = max(x(2:end)-x(1:end-1)); % delta x

%% Initial Conditions in physical Space
% Load Macroscopic Fugacity [z], Velocity[u] and Temperature[t] 
    [z0,u0,t0,p0,rho0,E0,tEnd,~] = SSBGK_IC1d(x,IC_case);
   
%% Microscopic Velocity Discretization (for DOM)
% That is to make coincide discrete microscopic velocities as quadrature abscisas 
% for integrating 'f', the probability distribution in other to recover our
% macroscopics properties.
switch quad

    case{1} % Newton Cotes Quadrature for D.O.M.
    V  = [-20,20];  % range: a to b
    nv = 200;       % nodes desired (may not the actual value)
    [v,w,k] = cotes_xw(V(1),V(2),nv,5); % Using Netwon Cotes Degree 5
    nv = length(v);     v = repmat(v,1,nx);     w = repmat(w,1,nx);
    J = ones(size(v)); % Jacobian = 1;
    
    case{2} % Gauss Hermite Quadrature for D.O.M.
    nv = 80;        % number of nodes selected (the actual value)
    [v,w] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    w = w.*exp(v.^2);   % weighting function for fixed GH quadrature
    v = repmat(v,1,nx);     w = repmat(w,1,nx);
    J = ones(size(v));  % Jacobian = 1;
    
    case{3} % Gauss Hermite Quadrature for Dynamic-D.O.M.
    nv = 15;         % number of nodes selected (the actual value)
    [v,ws] = GaussHermite(nv); % for integrating range: -inf to inf
    k = 1;          % quadrature constant.
    v = repmat(v,1,nx);     ws = repmat(ws,1,nx);
        
    otherwise
        error('Order must be between 1, 2 and 3');
end
    
%% Applying Discrete Ordinate Method on ICs:
[z,ux,t] = apply_DOM(z0,u0,t0,nv);  % Semi-classical IC
[p,~,~] = apply_DOM(p0,rho0,E0,nv); % Classical IC

%% Initial Equilibrium Distribution 'f0'
% Semiclassical Equilibrium distribuition function,
f0 = 1./(exp( (v-ux).^2 ./ t)./z + theta);

% Compute new v and J values if DDOM is used
if quad == 3 % if DDOM, update a, J and v.
    % Define dynamic reference,
    a = sqrt(t); v_star = (v-ux)./a;
    % function weights for Dynamic GH quadrature,
    w = ws.*exp(v_star.^2);
    % Jacobian for 1D case,
    J = sqrt(t);
end
    
% Plot IC of Distribution function, f, in Phase-Space:
if plot_figs == 1 %&& (quad == 1 || quad == 2)
   figure(1)
   surf(f0); grid on;
   xlabel('x - Spatial Domain'); 
   ylabel('v - Velocity Space');
   zlabel('f - Probability');
end

% Compute Initial Macroscopic Momemts:
n   = k*sum(J.* w .* f0);  % Density [n]
nux = k*sum(J.^2.* (v_star+ux) .* w .* f0);  % Macrospic moment [n*ux]
ne  = 0.5*k*sum(J.^3 .*( (J.*v_star+ux).^2 ).* w .* f0);  % Energy Density [n*E]
nE  = 0.5*k*sum(J.^3 .*( (v-ux./J).^2 ).* w .* f0);  % Internal energy [n*e]

% Initial macroscopic properties need to draw figure
E = nE;
rho = n;
    
%% Marching Scheme
% Initial time
time = 0;

tic
switch method
    case{1} % Upwind O(h)
        % Load IC
        f = f0;
        % Main loop
        while time < tEnd
            % Update time step
            dt = dx*CFL/max(max(abs(v))); % however v_max is always known
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
            
            % Compute equilibrium distribution for the current t_step
            f_eq = 1./(exp( (v-ux).^2 ./ t)./z + theta);
            
            % initialize variables
            u_eq = zeros(1,nx);
            u  = zeros(1,nx);
                            
            % (this part can be done in parallel!)
            for i = 1:nv  % for DDOM i = 1,2,3
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
                [F_left,F_right] = Upwindflux1d(u,v(i,:));
             
                % Compute next time step
                 u_next = u - dtdx*(F_right - F_left) - (dt/r_time)*(u-u_eq);
                 
                % BC
                u_next(1) = u_next(2);
                u_next(nx) = u_next(nx-1);

                % UPDATE info
                u = u_next;
                             
                % Going back to f
                f(i,:) = u(:);
            end
                                               
            % Compute macroscopic moments
            n   = k*sum(J.* w .* f);  % Density [n]
            nux = k*sum(J.^2.* (v_star+ux) .* w .* f);  % Macrospic moment [n*ux]
            ne  = 0.5*k*sum(J.^3.*( (J.*v_star+ux).^2 ).* w .* f);  % Energy Density [n*E]
            nE  = 0.5*k*sum(J.^3.*( (v-ux./J).^2 ).* w .* f);  % Internal energy [n*e]
            
            % Compute ux and t
            ux = nux./n;
            t = 4*ne./n;
                        
            % Apply DOM
            [~,ux,t] = apply_DOM(1,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Compute Dynamic quadrature requirements
            if quad == 3 % if DDOM, update a, J and v.
                % Define dynamic reference,
                a = sqrt(t); v_star = (v-ux)./a;
                % function weights for Dynamic GH quadrature,
                w = ws.*exp(v_star.^2);
                % Jacobian for 1D case,
                J = sqrt(t);
            end

            % Compute macroscopic moments (again)
            n   = k*sum(J.* w .* f);  % Density [n]
            nux = k*sum(J.^2.* (v_star+ux) .* w .* f);  % Macrospic moment [n*ux]
            ne  = 0.5*k*sum(J.^3.*( (J.*v_star+ux).^2 ).* w .* f);  % Energy Density [n*E]
            nE  = 0.5*k*sum(J.^3.*( (v-ux./J).^2 ).* w .* f);  % Internal energy [n*e]
            
            % Update macroscopic properties
            rho = n;
            ux = nux./n;
            p = ne;
            t = 4*ne./n;
            z = n./sqrt(pi.*t);
            E = ne + 0.5*ux.^2;
            
            % Apply DOM
            [z,ux,t] = apply_DOM(z,ux,t,nv); % Semi-classical variables
            %[p,rho,E] = apply_DOM(p,rho,E,nv); % Classical variables
            
            % Update figures
            drawnow
        end
        
    case{2} % TVD ... coming soon
        
    otherwise
        error('Order must be between 1 and 2');
end
toc
%% Results
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