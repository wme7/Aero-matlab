%% Discontinuous Galerkin (DG) 1D subroutine
% to solve the scalar advection equation:
%
% du/dt + dF(u)/dx = 0
%
% where $dF/du = a$
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.08.31

clear all; close all; clc;

%% Routine main parameters
    a = -0.5;   % Scalar velocity in x direction
   Ne = 10;     % Number of elements
    k = 4;      % order of polynomial used for approximation
  cfl = 0.6;    % Courant Number
t_end = 0.5;    % End time
    
%% Domain Parameters
  a_p = max(0,a);       % a^{+}
  a_m = min(0,a);       % a^{-}
    V = [0, 2];         % Domain range
   dv = (V(2)-V(1))/Ne; % Spatial step size
   dt = cfl*dv/abs(a);  % Time step size
 dtdx = dt/dv;          % precomputed to save some flops

%% Discretization of Domain
    vx = V(1):dv:V(2);  % Grid points
    t = 0: dt :t_end;   % time steps

%% Compute Elements
   Nx = Ne+1;          % Number of points
% Element to node connectivity for the one 1d case:
 EtoV = zeros(Ne,2);
for i = 1:Ne  EtoV(i,1) = i;  EtoV(i,2) = i+1; end
    
%% Initial Condition
  u_0 = ones(1,Nx);
  v_1 = ceil(2*Nx/5); % point v1 = 0.4 [-]
  v_2 = ceil(3*Nx/5); % point v2 = 0.6 [-]
  u_0(v_1:v_2) = 2;

%% Compute Gauss Legendre Points, Weigths, and Vandermonde matrix
r = lglnodes(k);
P  = Vandermonde1D(k, r); invP = inv(P);
Dr = Dmatrix1D(k, r, P);

%% Build points in elements
va = EtoV(:,1); vb = EtoV(:,2);
 x = ones(k+1,1)*vx(va) + 0.5*(r+1)*(vx(vb)-vx(va));

%% Lift

 
%% Load Initial Condition
u = u_0;

%% Main Loop

