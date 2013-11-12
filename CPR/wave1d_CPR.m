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
% Notes: Basic Scheme Implementation without RK integration method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.01; % CFL condition
tEnd = 1; % final time
K = 5; % degree of accuaracy
nE = 10; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=+1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' %Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 1],nE,'Legendre',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);

% Build Lagrange k-Polynomials
l = LagrangePolynomial(xgrid.solutionPoints);
L.lcoef = double(subs(l.lagrangePolynomial,-1));
L.rcoef = double(subs(l.lagrangePolynomial,1));
L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));

% IC
u0 = IC(x,3);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),0.9*min(min(u0)),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; it = 0;

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    % Plot u
    plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
    
    % compute fluxes in node coordinates
    f = flux(u);

    % Interpolate u and flux values at the boundaries of Ij
    switch xgrid.quadratureType
        case 'LGL'
            u_lbd = u(1,:);
            u_rbd = u(end,:);
            f_lbd = f(1,:);
            f_rbd = f(end,:);
        otherwise
            u_lbd = L.lcoef*u;
            u_rbd = L.rcoef*u;
            f_lbd = L.lcoef*f;
            f_rbd = L.rcoef*f;
    end
    % Build Numerical fluxes acroos faces
    u_pface = [u_lbd,0]; % + side 
    u_nface = [0,u_rbd]; % - side 

    % Apply Periodic BCs
    %u_nface(1) = u_nface(end); % left BD
    %u_pface(end) = u_pface(1); % right BD
    
    % Apply Neumann BCs
    u_nface(1) = u_pface(1); % left BD
    u_pface(end) = u_nface(end);%u_pface(end); % right BD

    % LF numerical flux
    alpha = max(max(abs(dflux(u)))); 
    nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
    nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

    % flux derivate
    df = L.dcoef*f;

    % Compute the derivate: F = f + gL*(nfluxL-f_bdL) + gR*(nfluxR-f_bdR)
    dF = df + dg.RR*(nfluxL - f_lbd) + dg.RL*(nfluxR - f_rbd);

    % next time info!
    u_next = u - dt*dF/J;
    
    % update info
    u = u_next;
    
    if rem(it,10) == 0
        drawnow;
    end
end