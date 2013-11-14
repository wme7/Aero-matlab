%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with CPR/FR
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: A flux reconstruction approach to high-order schemes including
% Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Simulation Parameters
fluxfun = 'nonlinear'; % select flux function
cfl = 0.005; % CFL condition
tEnd = 1.5; % final time
K = 3; % degree of accuaracy %example: K = 6 -> cfl 0.001
nE = 40; % number of elements
M = 10; % MODminmod parameter

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=-1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' %Burgers 
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 2*pi],nE,'LGL',K);
dx = xgrid.elementSize;     J = xgrid.Jacobian; 
x = xgrid.nodeCoordinates;  quad = xgrid.quadratureType;
w = xgrid.weights';     xc = xgrid.elementCenter;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);

% Build Lagrange k-Polynomials
l = LagrangePolynomial(xgrid.solutionPoints);
L.lcoef = double(subs(l.lagrangePolynomial,-1));
L.rcoef = double(subs(l.lagrangePolynomial,1));
L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));

% IC
ic = 2; u0 = IC(x,ic);
if ic==2
    load('burgersExact.mat');
    ue = burgersExact(2,:);
    xe = burgersExact(1,:);
end

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; it = 0;

% Using a 3-stage TVD Runge Kutta time integration
while t < tEnd
    uo = u;
    
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    % Limit solution
    [u,tcells] = limitSolution(u,xgrid,M);
    disp(tcells);
    
    % 1st stage
    dF = residual(u,L,dg,flux,dflux,quad);
    u = uo-dt*dF/J;

    % 2nd Stage
    dF = residual(u,L,dg,flux,dflux,quad); 
    u = 0.75*uo+0.25*(u-dt*dF/J);

    % 3rd stage
    dF = residual(u,L,dg,flux,dflux,quad); 
    u = (uo+2*(u-dt*dF/J))/3;
    
    % build cell averages
    u_bar = w*u/2;
    
    % Plot u
    subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange); grid on; 
    subplot(1,2,2); plot(xc,u_bar,'ro'); axis(plotrange); grid off; 

    %if rem(it,10)==0
        drawnow;
    %end
end    
%% Final Plot for IC 2
if ic==2
    subplot(1,2,1); plot(x,u,x,u0,'-+'); axis(plotrange); grid on;
    subplot(1,2,2); plot(xe,ue,'k-',xc,u_bar,'ro'); axis(plotrange);
end