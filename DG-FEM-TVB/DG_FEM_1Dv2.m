%% One-dimensional Conservation Law Solver 
% Computing Burgers solution using a DG-FEM routine with Polynomial
% truncation for the nonlinear terms. 
%
% for solving the following problem:
%
% u_t + f(u)_x = s(u)
%
% where f(u) = a*u  and  s(u) = u^2 
% Function residual will be defined as:
%
% u_t = Residue(u)
%
% Based on ideas of the following papers:
%
% 1. TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite Element
% Method for conservation laws II: General Framework. (1989)
% 2. Runge-Kutta Discontinuous Galerkin Method Using WENO Limiters. (2005)
%
% Coded by Manuel Diaz 2012.12.05

%% Clear Work Space 
clc; clear all; close all;

%% Simulation Parameters
k         = 3;      % Space order / Number of degress of freedom: 0 to k
np        = k+1;    % Number of points per Cell/Element
quadn     = 1;      % Element Integration quadrature to use: n={1,2,3,4}
rks       = 3;      % Time order of the Sheme/RK stages
flux_type = 3;      % {1}Roe, {2}Global LF, {3}LLF
a         = 0.5;    % scalar advection speed
cfl       = 0.3;    % Courant Number
tEnd      = pi/15;  % Final Time for computation
nx        = 80;     % Number of Cells/Elements
M         = 0.01;   % TVB constant M
IC_case   = 3;      % 4 cases are available
plot_figs = 1;      % {1}Plot figures, {0}Do NOT plot figures

%% Define Grid Cells (Global) Points
% Building an uniform grid of cells/Elements: 
% Grid points are defined as the center points of
% x_{j} and x_{j+1}. That is x_{j+1/2}.
x_left = 0; x_right = 2; dx = (x_right-x_left)/nx;
x_temp    = x_left: dx :x_right;
x         = repmat(x_temp,np,1); % Cells

%% Build Cells/Elements (Local) Points
switch quadn
    case{1} % Scaled Guass-Lobatto quadrature weights and abscissas 
        [xi,w]   = sGaussLobatto(np); % scaled to [-1/2, 1/2]
        xj       = repmat( xi*dx ,1,nx+1);
        w        = repmat(   w   ,1,nx+1);
        x        = x + xj; % Create poinst in every cells
        
    case{2} % Guass-Lobatto quadrature weights and abscissas 
        [xi,w,V] = GaussLobatto(np); % for the interval [-1,1]
        xj       = repmat( xi*dx ,1,nx+1);
        w        = repmat(   w   ,1,nx+1);
        x        = x + xj; % Create poinst in every cells
        
    case{3} % Guass-Legendre quadrature weights and abscissas
        [xi,w] = GaussLegendre(np); % for the interval [-1,1]
        xj       = repmat( xi*dx ,1,nx+1);
        w        = repmat(   w   ,1,nx+1);
        x        = x + xj; % Create poinst in every cells
                
    case{4} % Guass-Legendre quadrature weights and abscissas 
        [xi,w,V] = GaussRadau(np); % for the interval [-1,1]
        xj       = repmat( xi*dx ,1,nx+1);
        w        = repmat(   w   ,1,nx+1);
        x        = x + xj; % Create poinst in every cells
        
    otherwise
        error('quadrature not available')
end
%% Computing Legendre Polynomials 
% We use them as basis for expanding our solution into degress of freedom



%% flux function
flux = @(c,u) c*u;

%% source term fucntion
Sterm = @(u) u.^2;

%% Load Initial Condition, u(x,0) = u0
u0 = u_zero(x,IC_case);
f0 = flux(a,u0);
s0 = Sterm(u0);

% Plot u0
%if plot_figs == 1; plot(x,u0); axis([0,2,0,1.2]); end;
if plot_figs == 1; 
    subplot(1,3,1); plot(x,u0); title('u0'); axis tight;
    subplot(1,3,2); plot(x,f0); title('f0'); axis tight;
    subplot(1,3,3); plot(x,s0); title('s0'); axis tight; 
end;

% Define coeficients matrix
al = diag([1 12 180 2800 44100 698544 11099088 176679360]);

% Transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
P = sLegMat(k,xi);   % Legendre Matrix

ut = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx+1
        ut(i,j) = al(i,i).*sum(w(:,j).*u0(:,j).*P(:,i));
    end
end

% Transform f(u) to degress of freedom f(t)_{l,i} for each i-Cell/Element
ft = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx+1
        ft(i,j) = al(i,i).*sum(w(:,j).*f0(:,j).*P(:,i));
    end
end

% Transform s(u) to degress of freedom s(t)_{l,i} for each i-Cell/Element
st = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    i = l+1;            % Dummy index
    for j = 1:nx+1
        st(i,j) = al(i,i).*sum(w(:,j).*s0(:,j).*P(:,i));
    end
end

%% Elements boundary values
%fp = f( 1,:); % flux to the left  (f( 1,i)^{+})
%fn = f(np,:); % flux to the right (f(np,1)^{-})

%% Compute Residue
% STEP 1
%-------
% Cell Integral volume by Gauss-Lobatto quadrature.
% NOTE:
%  1. This is the long way: element by element.
%  2. This is very slow because of the two 'for' loops.
v_residue = zeros(np,nx+1);
for l = 0:k
    i = l+1;    % Dummy index
    for j = 1:nx
        v_residue(i,j) = al(i,i)*sum(w(:,j).*flux(a,u0(i,j)).*dsLegendreP(l,xi));
    end
end

% Cell Integral volume by Prof. Tim Warburton method.
% NOTE:
% 1. this is a faster and very mechanic way to compute this term.


% STEP 2
%-------
% Cell boundary contribution
b_residue = zeros(np,nx+1);
up = u0( 1,:); % state at left boundary (u^{+} = u( 1,i))
un = u0(np,:); % state at right boundary (u^{-} = u(np,i))


%% Transform degress of freedom u(t)_{l,i} back to space-time values u(x,t)
u = (ut'*P')';
f = a*u;
s = u.^2;
figure
if plot_figs == 1; 
    subplot(1,3,1); plot(x,u0); title('u'); axis tight;
    subplot(1,3,2); plot(x,f0); title('f'); axis tight;
    subplot(1,3,3); plot(x,s0); title('s'); axis tight; 
end;