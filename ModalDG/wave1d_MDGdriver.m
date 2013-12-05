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
% Notes: Basic Scheme Implementation with SSP-RK33 intergration scheeme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.01; % CFL condition
tEnd = 2*pi; % final time
K = 7; % degree of accuaracy
nE = 10; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=-1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' % Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 2*pi],nE,'LGL',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;
w = xgrid.weights';	xc = xgrid.elementCenter;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vandermonde;
lR = tool.legRightEnd; lL = tool.legLeftEnd;
D = tool.MassMatrix*tool.CoefDiffMatrix;
invM = tool.invMassMatrix;

% IC
u0 = IC(x,2);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; ut = V\u; it = 0;

% Using a 3-stage TVD Runge Kutta time integration
while t < tEnd
    uo = ut;
    
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    % 1st stage
    dF = residual(u,ut,flux,dflux,lR,lL,D,V,invM);
    ut = uo-dt*dF/J;

    % 2nd Stage
    dF = residual(u,ut,flux,dflux,lR,lL,D,V,invM); 
    ut = 0.75*uo+0.25*(ut-dt*dF/J);

    % 3rd stage
    dF = residual(u,ut,flux,dflux,lR,lL,D,V,invM); 
    ut = (uo+2*(ut-dt*dF/J))/3;
    
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
subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange);
title('Modal DG','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);
subplot(1,2,2); plot(x,u0,'k-',xc,u_bar,'ro'); axis(plotrange);
title('Cell Averages','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);