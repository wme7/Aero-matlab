%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with CPR/FR
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Simulation Parameters
fluxfun = 'nonlinear'; % select flux function
cfl = 0.01; % CFL condition
tEnd = 2; % final time
K = 5; % degree of accuaracy
nE = 10; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=-2; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' %Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh
xgrid = mesh1d([0 2],nE,'LGL',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dRR = RR.eval_dP(xgrid.solutionPoints); dRL = -flipud(dRR);

% Build Lagrange k-Polynomials
L = LagrangePolynomial(xgrid.solutionPoints);
L_lcoef = double(subs(L.lagrangePolynomial,-1));
L_rcoef = double(subs(L.lagrangePolynomial,1));
dL_coef = double(subs(L.dlagrangePolynomial,xgrid.solutionPoints));

% IC
u0 = IC(x,2);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),0.9*min(min(u0)),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; it = 0;

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u))));
    t = t+dt;
    
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
            u_lbd = L_lcoef*u;
            u_rbd = L_rcoef*u;
            f_lbd = L_lcoef*f;
            f_rbd = L_rcoef*f;
    end
    % Build Numerical fluxes at faces
    u_faces(1,:) = [u_lbd,0];
    u_faces(2,:) = [0,u_rbd];
    
    %u_pface = [u_lbd,0];
    %u_nface = [0,u_rbd];

    % Apply Periodic BCs
    u_faces(2,1) = u_faces(2,end); % left BD
    u_faces(1,end) = u_faces(1,1); % right BD
    %u_nface(1) = u_nface(end); % left BD
    %u_pface(end) = u_pface(1); % right BD

    % evaluate LF
    alpha = max(max(abs(dflux(u)))); 
    LF = @(l,r) 0.5*(flux(l)+flux(r)-alpha*(r-l));
    nflux = LF(u_faces(2,:),u_faces(1,:));
    %nflux = 0.5*(flux(u_nface)-flux(u_pface)-alpha*(u_pface-u_nface));
    %nflux = 0.5*(flux(u_pface)-flux(u_nface)-alpha*(u_nface-u_pface));
    nfluxL = nflux(1:end-1);
    nfluxR = nflux(2:end);

    % flux derivate
    df = dL_coef*f;

    % Compute the derivate: F = f + g*(nfluxL - f_bdL) + g*(nfluxR - f_bdR)
    dF = df + dRR*(nfluxL - f_lbd) + dRL*(nfluxR - f_rbd);

    % next time info!
    u_next = u - dt*dF/J;
    
    % update info
    u = u_next;
    
    if rem(it,10) == 0
        drawnow;
    end
end