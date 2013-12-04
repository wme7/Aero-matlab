%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with Nodal DG
%
%               du/dt + df/dx = 0,  for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2012.12.05
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: J.S. Hestaven & T. Warburton; Nodal Discontinuous Galerkin Methods,
% Algorithms, Analysis and Applications. Springer 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Scheme Implementation without RK intergration scheeme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; %clc;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.02; % CFL condition
tEnd = 6*pi; % final time
K = 5; % degree of accuaracy
nE = 10; % number of elements

%% PreProcess
% Define our Flux function
switch fluxfun
    case 'linear'
        a=1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    case 'nonlinear' % Burgers
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

% Build 1d mesh (LGL or ChebyshebMod only!)
xgrid = mesh1d([0 2*pi],nE,'LGL',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vandermonde2; 
invM = tool.nodalInvMassMatrix;
Dr = tool.nodalCoefDiffMatrix;

% Build Lift Operator
Emat = zeros(K+1,2); % array of element's shape function
Emat(1,1)=1; Emat(K+1,2)=1;
Lift = V*(V'*Emat);

% IC
u0 = IC(x,1);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t=0; u=u0; it=0;

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;

    % iteration counter
    it = it+1; 

    % compute fluxes in node coordinates
    f = flux(u); 

    % u and flux values at the boundaries of Ij
    u_lbd = u(1,:);
    u_rbd = u(end,:);
    f_lbd = f(1,:);
    f_rbd = f(end,:);

    % Build Numerical fluxes across faces
    u_pface = [u_lbd,0]; % + side 
    u_nface = [0,u_rbd]; % - side 

    % Apply Periodic BCs
    %u_nface(1) = u_nface(end); % left BD
    %u_pface(end) = u_pface(1); % right BD
    
    % Apply Neumann BCs
    u_nface(1) = u_pface(1); % left BD
    u_pface(end) = u_nface(end);% u_pface(end); % right BD
    
    % LF numerical flux
    alpha = max(max(abs(dflux(u)))); 
    nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
    nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);
    
    % Compute the derivate: F = f + lL*(nfluxL-f_bdL) - lR*(nfluxR-f_bdR)
    dF = -Dr*f + Lift*[(nfluxL-f_lbd);(f_rbd-nfluxR)];
    
    % next time info!
    ut_next = u + dt*dF/J;

    % UPDATE info
    u = ut_next;

    % Plot u
    %if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
        xlabel('x'); ylabel('u'); title('DG-FEM')
        drawnow;
    %end
end