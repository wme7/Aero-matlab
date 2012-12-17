%% One-dimensional Conservation Law Solver 
% Computing Burgers solution using a DG-FEM routine with Polynomial
% truncation for the nonlinear terms. 
%
% for solving the following problem:
%
% u_t + f(u)_x = s(u)
%
% Based on ideas of the following papers:
%
% 1. TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite Element
% Method for conservation laws II: General Framework. (1989)
% 2. Runge-Kutta Discontinuous Galerkin Method Using WENO Limiters. (2005)
%
% Coded by Manuel Diaz 2012.12.05

%% Clear Work Space 
clear all; close all; %clc;

%% Simulation Parameters
k         = 3;      % Space order / Number of degress of freedom: 0 to k
np        = k+1;    % Number of points per Cell/Element
quadn     = 3;      % Element Integration quadrature to use: n={1,2,3,4}
rks       = 3;      % Time order of the Sheme/RK stages
flux_type = 1;      % {1}Roe, {2}Global LF, {3}LLF
cfl       = 0.3;    % Courant Number
tEnd      = pi/15;  % Final Time for computation
nx        = 10;     % Number of Cells/Elements
MM        = 0.01;   % TVB constant M
IC_case   = 3;      % 4 cases are available
plot_figs = 1;      % {1}Plot figures, {0}Do NOT plot figures

%% Define Grid Cells (Global) Points
% Building an uniform grid of cells/Elements: 
% Grid points are defined as the center points of
% x_{j} and x_{j+1}. That is x_{j+1/2}.
x_left = 0; x_right = 1; dx = (x_right-x_left)/nx;
x_temp    = x_left+dx/2: dx :x_right-dx/2;
x         = repmat(x_temp,np,1); % Cells

%% Build Cells/Elements (Local) Points
switch quadn
    case{1} % Scaled Guass-Lobatto quadrature weights and abscissas 
        [xi,w]   = sGaussLobatto(np); % scaled to [-1/2, 1/2]
        xj       = repmat( xi*dx ,1,nx);
        w        = repmat(   w   ,1,nx);
        x        = x + xj; % Create poinst in every cells
        
    case{2} % Guass-Lobatto quadrature weights and abscissas 
        [xi,w,V] = GaussLobatto(np); % for the interval [-1,1]
        xj       = repmat( xi*dx/2 ,1,nx);
        w        = repmat(   w   ,1,nx);
        x        = x + xj; % Create poinst in every cells
        
    case{3} % Guass-Legendre quadrature weights and abscissas
        [xi,w] = GaussLegendre(np); % for the interval [-1,1]
        xj       = repmat( xi*dx/2 ,1,nx);
        w        = repmat(   w   ,1,nx);
        x        = x + xj; % Create poinst in every cells
                
    case{4} % Guass-Legendre quadrature weights and abscissas 
        [xi,w,V] = GaussRadau(np); % for the interval [-1,1]
        xj       = repmat( xi*dx/2 ,1,nx);
        w        = repmat(   w   ,1,nx);
        x        = x + xj; % Create poinst in every cells
        
    otherwise
        error('quadrature not available')
end

%% Load Initial Condition, u(x,0) = u0
u0 = u_zero(x,IC_case);

%% Math Objects
M = zeros(np,np);
D = zeros(np,np);
C = zeros( 1,np);
for l = 0:k             % For all degress of freedom
    i = l+1;            % Dummy index
    for j=1:np          % For all local points
        if j>i && rem(j-i,2)==1
            D(i,j)=2;   % D or diffentiated Matrix
        end
    end
    M(i,i) = (2*l+1)/dx; % Inverser of Mass Matrix
    C(1,i) = (-1)^l;    % Legendre polynomials of deg 'l' evaluated at x=1
end

%% Plot u0
if plot_figs == 1; plot(x,u0); axis auto; grid on; end;

% Pre-allocate math objects
P = LegMat(k,xi);     % Legendre Vandermonde Matrix
ut = zeros(np,nx);    % u degress of freedom
St = zeros(np,nx);    % source degress of freedom
Ft = zeros(np,nx);    % flux degress of freedom

%% Computing the Evolution of the residue 'L', du/dt = L(u) 
% Load Initial condition
u = u0;
% Set Initial time step
time = 0; % time
n    = 0; % counter
%while time < tEnd
    dt = cfl*dx/max(max(abs(u)));   % dt, time step
    time = time + dt;               % iteration time
    n    = n + 1;                   % update counter

% transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
%ut = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx
        ut(i,j) = (2*l+1)/2*sum(w(:,j).*u(:,j).*P(:,i));
    end
end
    
% Source Term S(u) = (u)^2
S  = u.^2;
%St = zeros(np,nx); 
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx
        St(i,j) = (2*l+1)/2.*sum(w(:,j).*S(:,j).*P(:,i));
    end
end

% Flux Term f(u) = (u)^2/2 thus f'(u) = u
F  =  1*u; % u.^2/2;
Ft = zeros(np,nx);
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx
        Ft(i,j) = (2*l+1)/2.*sum(w(:,j).*F(:,j).*P(:,i));
    end
end

% Volume Term contribution
volume = D'*Ft;

% Fluxes contribution


% We use simple Upwinding Scheme
ut_next = ut;

% Transfrom back from degress of freedom
u = (ut'*P')';
%end % while-end

%% Transform degress of freedom u(t)_{l,i} back to space-time values u(x,t)
%u = (ut'*P')';
figure
if plot_figs == 1; plot(x,u); axis tight; grid on; end;
