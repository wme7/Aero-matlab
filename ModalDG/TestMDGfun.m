function [L1_nodes,degree_aveg,Linf_nodes,degree_min] = TestMDGfun(fluxfun,cfl,tEnd,Ic,K,nE,quad)
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
%clear all; close all; %clc;

%% Parameters
% fluxfun = 'linear'; % select flux function
% cfl = 0.02; % CFL condition
% tEnd = 2*pi; % final time
% K = 3; % degree of accuaracy
% nE = 20; % number of elements
%quad = 'LGL';

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
xgrid = mesh1d([0 2*pi],nE,quad,K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;
%w = xgrid.weights';	xc = xgrid.elementCenter;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vandermonde;
lR = tool.legRightEnd; lL = tool.legLeftEnd;
D = tool.MassMatrix*tool.CoefDiffMatrix;
invM = tool.invMassMatrix;

% IC
u0 = IC(x,Ic);

% Set plot range
%plotrange = [xgrid.range(1),xgrid.range(2),...
%    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; ut = V\u; it = 0;

% Using a 3-stage TVD Runge Kutta time integration
% Fixed time evolution
dt = cfl*dx/abs(a); 
time = t:dt:tEnd;
%while t < tEnd
for s = time
    uo = ut;
    
    % update time
    %dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
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
    %u_bar = w*u/2;
    
    % Plot u
    %subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange); grid off; 
    %subplot(1,2,2); plot(xc,u_bar,'ro'); axis(plotrange); grid off; 
    
    %if rem(it,10) == 0
    %    drawnow;
    %end
end

%% Build Cell averages for every E_j
%u0_avg = w*u0;
%u_avg = w*u;

%% Compute Norms
err = u0-u;

% L1 Error
%L1_avgs = sum(abs(u0_avg-u_avg));
L1_nodes = sum(sum(abs(err)));
degree_aveg = abs(log10(L1_nodes/xgrid.nNodes));

% L\infty Error
%Linf_avgs = max(abs(u0_avg-u_avg));
Linf_nodes = max(max(abs(err)));
degree_min = abs(log10(Linf_nodes));