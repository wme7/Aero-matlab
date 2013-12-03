%% Driver script for solving the 1D advection equations
Globals1D;

%% Order of polymomials used for approximation 
N = 3;

%% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,1,10);

%% Initialize solver and construct grid and metric
StartUp1D;

%% Set initial conditions
u = sin(2*pi*x);

%% Solve Problem
FinalTime = 8*pi;
[u] = Advec1D(u,FinalTime);

plot(x,u);