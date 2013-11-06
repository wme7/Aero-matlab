%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Solving 1-D scalar wave equation with CPR/FR
%
%                       du/dt + a dudx = 0
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;

%% Simulation Parameters
a = -1; % advection speed
cfl = 0.1; % CFL condition
tEnd = 0.2; % final time
K = 4; % degree of accuaracy
nE = 8; % number of elements

%% Prepocessing 
% Define our Flux function
flux = @(w) a*w; %w.^2/2; 

% Derivate of the flux function
dflux = @(w) a*ones(size(w)); %w;

% Build 1d mesh
grid = mesh1d([0 1],nE,'Legendre',K);
dx = grid.elementSize; J = grid.Jacobian; x = grid.nodeCoordinates;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dRR = RR.eval_dP(grid.solutionPoints); dRL = -flipud(dRR);

% Build Lagrange k-Polynomials
L = LagrangePolynomial(grid.solutionPoints);

% IC
u0 = IC(x,1);

%% Solver Loop
% load IC
u = u0;

% time step
dt = cfl*dx/abs(a);
dtdx = dt/dx;

for time = 0:dt:tEnd

    % Plot u
    plot(x,u,'-o')

    % compute fluxes in node coordinates
    f = flux(u);

    % Interpolate u at the boundaries of [-1,1]
    u_bd = zeros(2,grid.nElements);
    u_bd(1,:) = double(subs(L.lagrangePolynomial,-1))*u;
    u_bd(2,:) = double(subs(L.lagrangePolynomial,1))*u;

    u_faces = zeros(2,grid.nFaces);
    u_faces(1,:) = [u_bd(1,:),0];
    u_faces(2,:) = [0,u_bd(2,:)];

    % Periodic BCs
    u_faces(2,1) = u_faces(2,end); % left BD
    u_faces(1,end) = u_faces(1,1); % right BD

    % Numerical Flux: LF
    alpha = max(max(abs(dflux(u)))); 
    LF = @(a,b) 0.5*(flux(a)+flux(b)-alpha*(b-a));
    nflux = LF(u_faces(1,:),u_faces(2,:));
    nfluxL = nflux(1:end-1);
    nfluxR = nflux(2:end);

    % interpolate flux at the boundaries of [-1,1]
    f_bdL = double(subs(L.lagrangePolynomial,-1))*f;
    f_bdR = double(subs(L.lagrangePolynomial,1))*f;

    % flux derivate
    df = double(subs(L.dlagrangePolynomial,grid.solutionPoints))*f;

    % Compute the derivate: F = f + g*(nfluxL - f_bdL) + g*(nfluxR - f_bdR)
    dF = df + dRR*(nfluxL - f_bdL) + dRL*(nfluxR - f_bdR);

    % next time info!
    u_next = u - dt*J*dF;
    
    % update info
    u = u_next;

    % update plot
    drawnow

end