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
% Notes: Basic Scheme Implementation without RK integration method and no
% 'for' loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; %clc;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.02; % CFL condition
tEnd = 2*pi; % final time
K = 5; % degree of accuaracy
nE = 10; % number of elements

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

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vadermonde; lR = tool.legRightEnd; lL = tool.legLeftEnd;
M = tool.MassMatrix; invM = tool.invMassMatrix; Dr = tool.CoefDiffMatrix;
D = M*Dr;

% IC
u0 = IC(x,1);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; ut = V\u; it = 0;

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
    
    % compute fluxes in node coordinates
    f = flux(u); ft = V\f;

    % Interpolate u and flux values at the boundaries of Ij
    u_lbd = lL*ut;
    u_rbd = lR*ut;
    
    % Build Numerical fluxes across faces
    u_pface = [u_lbd,0]; % + side 
    u_nface = [0,u_rbd]; % - side 

    % Apply Periodic BCs
    u_nface(1) = u_nface(end); % left BD
    u_pface(end) = u_pface(1); % right BD

    % Apply Neumann BCs
    %u_nface(1) = u_pface(1); % left BD
    %u_pface(end) = u_nface(end);% u_pface(end); % right BD

    % LF numerical flux
    alpha = max(max(abs(dflux(u)))); 
    nflux = 0.5*(flux(u_nface)+flux(u_pface)-alpha*(u_pface-u_nface));
    nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

    % Compute the derivate: F = f - lR*(nfluxR) + lL*(nfluxL)
    dF = invM*(D'*ft - lR'*nfluxR + lL'*nfluxL);

    % next time info!
    ut_next = ut + dt*dF/J;
    
    % UPDATE info
    ut = ut_next;
       
    % Plot u
    %if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
        xlabel('x'); ylabel('u'); title('DG-FEM')
        drawnow;
    %end
    
    % Transform legendre coefs into nodal values.
    u = V*ut;
end