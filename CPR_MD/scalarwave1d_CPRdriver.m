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
cfl = 0.05; % CFL condition
tEnd = 1; % final time
K = 5; % degree of accuaracy
nE = 8; % number of elements

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
grid = mesh1d([0 2],nE,'LGL',K);
dx = grid.elementSize; J = grid.Jacobian; x = grid.nodeCoordinates;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dRR = RR.eval_dP(grid.solutionPoints); dRL = -flipud(dRR);

% Build Lagrange k-Polynomials
L = LagrangePolynomial(grid.solutionPoints);

% IC
u0 = IC(x,2);

%% Solver Loop

% Set initial time & load IC
t = 0; u = u0; 

while t < tEnd
    % update time
    dt = cfl*dx/max(max(abs(dflux(u))));
    t = t + dt;
    
    % Plot u
    plot(x,u,'-o')

    % compute fluxes in node coordinates
    f = flux(u);

    % Interpolate u at the boundaries of [-1,1]
    switch grid.quadratureType
        case 'LGL'
            u_bd = [u(1,:);u(end,:)];
        otherwise
            u_bd = zeros(2,grid.nElements);
            u_bd(1,:) = double(subs(L.lagrangePolynomial,-1))*u;
            u_bd(2,:) = double(subs(L.lagrangePolynomial,1))*u;
    end
 
    u_faces = zeros(2,grid.nFaces);
    u_faces(1,:) = [u_bd(1,:),0];
    u_faces(2,:) = [0,u_bd(2,:)];

    % Periodic BCs
    u_faces(2,1) = u_faces(2,end); % left BD
    u_faces(1,end) = u_faces(1,1); % right BD

    % Numerical Flux: LF
    alpha = max(max(abs(dflux(u)))); 
    LF = @(a,b) 0.5*(flux(a)+flux(b)-alpha*(b-a));
    nflux = LF(u_faces(2,:),u_faces(1,:));
    nfluxL = nflux(1:end-1);
    nfluxR = nflux(2:end);

    % interpolate flux at the boundaries of [-1,1]
    switch grid.quadratureType
        case 'LGL'
            f_bdL = f(1,:);
            f_bdR = f(end,:);
        otherwise
            f_bdL = double(subs(L.lagrangePolynomial,-1))*f;
            f_bdR = double(subs(L.lagrangePolynomial,1))*f;
    end

    % flux derivate
    df = double(subs(L.dlagrangePolynomial,grid.solutionPoints))*f;

    % Compute the derivate: F = f + g*(nfluxL - f_bdL) + g*(nfluxR - f_bdR)
    dF = df + dRR*(nfluxL - f_bdL) + dRL*(nfluxR - f_bdR);

    % next time info!
    u_next = u - dt*dF/J;
    
    % update info
    u = u_next;

    % update plot
    drawnow

end