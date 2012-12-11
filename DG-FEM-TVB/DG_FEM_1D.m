%% One-dimensional Conservation Laws Solver 
% Computing Burgers solution using a DG-FEM routine with Polynomial
% truncation for the nonlinear terms. 
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
mt        = 3;      % Time order of the Sheme
flux_type = 3;      % {1}Roe, {2}Global LF, {3}LLF
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
% Compute Guass-Lobatto quadrature weights and abscissas for every cell:
[xi,w]   = GaussLobatto(np); % this are scaled Gauss Lobatto poinst and weights
xj       = repmat( xi*dx ,1,nx+1);
w        = repmat(   w   ,1,nx+1);
x        = x + xj; % Create poinst in every cells

%% Load Initial Condition, u(x,0) = u0
u0 = u_zero(x,IC_case);

% Plot u0
%if plot_figs == 1; plot(x,u0); axis([0,2,0,1.2]); end;
if plot_figs == 1; plot(x,u0); axis('tight'); end;

% Define coeficients matrix
a = diag([1 12 180 2800 44100 698544 11099088 176679360]);

% transform u(x,t) to degress of freedom u(t)_{l,i} for each i-Cell/Element
ut = zeros(np,nx+1);
for l = 0:k             % for all degress of freedom
    for j = 1:nx
    i = l+1;            % Dummy index
    P = LegMat(k,xi);   % Legendre Matrix
    ut(i,j) = a(i,i).*sum(w(:,j).*u0(:,j).*P(:,i));
    end
end

%% transform degress of freedom u(t)_{l,i} back to space-time values u(x,t)
u = (ut'*P')';
figure
plot(x,u); axis('tight');