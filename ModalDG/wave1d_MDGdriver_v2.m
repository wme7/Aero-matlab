%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with Modal DG
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2012.12.05
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite
% Element Method for conservation laws II: General Framework. (1989)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Scheme Implementation with SSP-RK45 intergration scheeme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
fluxfun = 'nonlinear'; % select flux function
cfl = 0.02; % CFL condition
tEnd = 1.5; % final time
K = 5; % degree of accuaracy
nE = 20; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=+1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' % Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 2*pi],nE,'Radau',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;
w = xgrid.weights';	xc = xgrid.elementCenter;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vadermonde;

% IC
ic = 2; u0 = IC(x,ic);
switch ic
    case 2
    load('burgersExact.mat');
    ue = burgersExact(2,:);
    xe = burgersExact(1,:);
    case 3
    load('burgersExact2.mat');
    ue = burgersExact(2,:);
    xe = burgersExact(1,:);
    otherwise
        % do nothing
end

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Set initial time & load IC
t = 0; u = u0; ut = V\u; it = 0;

% Using a 4rd Order 5-stage SSPRK time integration
res_ut = zeros(K+1,nE); % Runge-Kutta residual storage

% Using a 3-stage TVD Runge Kutta time integration
while t < tEnd
    uo = ut;
    
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        dF = residual(u,ut,flux,dflux,tool);
        res_ut = rk4a(RKs)*res_ut + dt*dF/J;
        ut = ut - rk4b(RKs)*res_ut;
    end
    
    % Transform legendre coefs into nodal values.
    u = V*ut;
    
    % build cell averages
    u_bar = w*u/2;
    
    % Plot u
    subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange); grid off; 
    subplot(1,2,2); plot(xc,u_bar,'ro'); axis(plotrange); grid off; 
    
    %if rem(it,10) == 0
        drawnow;
    %end
end
%% Final Plot for IC 2
if ic==2 || ic==3
    subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange);
    title('Modal DG','interpreter','latex','FontSize',18);
    xlabel('$\it{x}$','interpreter','latex','FontSize',14);
    ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);
    subplot(1,2,2); plot(xe,ue,'k-',xc,u_bar,'ro'); axis(plotrange);
    title('Cell Averages','interpreter','latex','FontSize',18);
    xlabel('$\it{x}$','interpreter','latex','FontSize',14);
    ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);
end