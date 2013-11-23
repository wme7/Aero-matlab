%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with Modal DG
%
%               du/dt + df/dx = S(u), for x \in [a,b]
%                 where f = f(u): linear/nonlinear
%
%              coded by Manuel Diaz, NTU, 2012.12.05
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite
% Element Method for conservation laws II: General Framework. (1989)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Scheme Implementation without RK integration method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; %clc;

%% Parameters
fluxfun = 'linear'; % select flux function
cfl = 0.02; % CFL condition
tEnd = 2*pi; % final time
K = 3; % degree of accuaracy
nE = 20; % number of elements
includeS = 0; %{0}no sourcer term, {1}add source term
flux_type = 4; %{1}Roe, {2}LF, {3}LLF, {4}Upwind

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

% Source term function
switch includeS
    case{0} % no source term
        S = @(w) zeros(size(w));
    case{1} % with source term
        % example source term
        S = @(w) w.^2; 
end

% Build 1d mesh
xgrid = mesh1d([0 2*pi],nE,'LGL',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; x = xgrid.nodeCoordinates;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vadermonde; lR = tool.legRightEnd; lL = tool.legLeftEnd;
M = tool.MassMatrix; invM = tool.invMassMatrix; Dr = tool.CoefDiffMatrix;
D = M*Dr;

% IC
u0 = IC(x,2);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*u0)),min(min(1.1*u0))),1.1*max(max(u0))];

%% Solver Loop

% Load Initial conditions
u = u0; ut = V\u; t = 0; it = 0;

residue = zeros(K+1,nE);
while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u)))); t = t+dt;
    
    % iteration counter
    it = it+1; 
           
    % Evaluate/update f and s terms on domain
    f = flux(u); ft = V\f;
    s = S(u); st = V\s;

    % Compute fluxes - Inner cells boundaries
    up = lL*ut; % u_{i-1/2}^(+) -> Left u
    un = lR*ut; % u_{i+1/2}^(-) -> Right u
    uc = [un(1:end-1);up(2:end)];
    h  = tool.DGnflux(flux,dflux,uc,flux_type); % Evaluate fluxes
    
    % Compute fluxes for periodic BC - boundary cells
    ub = [un(end);up(1)];
    hb = tool.DGnflux(flux,dflux,ub,flux_type);
    
    % Just for comparison with mathematica
    %totalh = [hb,h,hb];
    
    % Step-by-step procedure:
    % Periodic BC @ cell#1
         residue(:,1) = ( D'*ft(:,1) - ...          % Volume term
                        h(1)*lR' + hb*lL'  ...       % flux terms
                        ).*diag(invM) + st(:,1);      % Source term
    
    % Compute residue function:
    for i = 2:nE-1
        residue(:,i) = ( D'*ft(:,i) - ...           % Volume term
                        h(i)*lR' + h(i-1)*lL'  ...   % flux terms
                        ).*diag(invM) + st(:,i);      % Source term
    end
    
    % Periodic BC @ cell#nx
        residue(:,nE) = ( D'*ft(:,nE) - ...         % Volume term
                        hb*lR' + h(nE-1)*lL'  ...    % flux terms
                        ).*diag(invM) + st(:,nE);     % Source term
    
    % Compute next time step
    ut_next = ut + dt*residue/J;

    % UPDATE info
    ut = ut_next;
    
   % Plot 
    %if rem(it,10) == 0
        plot(x,u0,'-x',x,u,'-'); axis(plotrange); grid on; 
        xlabel('x'); ylabel('u'); title('DG-FEM')
        drawnow;
    %end

    % Transform legendre coefs into nodal values.
    u = V*ut;
end