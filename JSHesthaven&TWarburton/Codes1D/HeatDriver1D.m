% Driver script for solving the 1D advection equations with variable coefficient

Globals1D;

% Polynomial order used for approximation 
N = 4;

% Read in Mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
%u0 = sin(x);

u0 = 0.1*ones(size(x));
x1 = -1/6; x2 = +1/6;
rhs = find(x>=x1 & x<=x2);
u0(rhs) = 1;

% Solve Problem
FinalTime = 0.05;
[u,time] = Heat1D(u0,FinalTime);
