%% One-dimensional Conservation Laws Solver 
% Computing Burgers solution using a DG-FEM routine with Polynomial
% truncation for the nonlinear terms. 
%
% Based on:
%
% 1. TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite Element
% Method for conservation laws II: General Framework. (1989)
% 2. Runge-Kutta Discontinuous Galerkin Method Using WENO Limiters. (2005)

%% Clear Work Space 
clc; clear all; close all;

%% Simulation Parameters
k     = 3;      % Space order for the Scheme
np    = k+1;    % Number of points per Cell/Element
mt    = 3;      % Time order of the Sheme
kflux = 3;      % {1}Roe, {2}Global LF, {3}LLF
cfl   = 0.3;    % Courant Number
tEnd  = pi/15;  % Final Time for computation
nx    = 80;     % Number of Cells/Elements
M     = 0.01;   % TVB constant M
show  = 1;      % {1}Plot figures, {2}Do NOT plot figures

%% Define Grid Cells (Global) Points
% Building an uniform grid of cells/Elements: 
% Grid points are defined as the center points of
% x_{j} and x_{j+1}. That is x_{j+1/2}.
x_left = 0; x_right = 2; dx = (x_right-x_left)/nx;
temp = (x_left+dx/2):dx:(x_right-dx/2);
x = repmat(temp,np,1); % Cells

%% Build Cells/Elements (Local) Points
% Compute Guass-Lobatto quadrature weights and abscissas for every cell:
[Cx,Cw] = GaussLobatto(np);
Cx = repmat( Cx*dx ,1,nx);
Cw = repmat(  Cw   ,1,nx);
x = x + Cx; % Create poinst in every cells

%% Load Initial Condition, u(x,0) 
u0 = 1/2+sin(pi*x);
% Plot u(x,0)
if show == 1; plot(x,u0); axis('tight'); end;

% Legendre Quadrature
u = zeros(np,nx);
for i = 0
    u(i+1,:) = sum(Cw.*u0.*LegendreP(i,x));
end

