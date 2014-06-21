%% Driver script for solving the 1D advection equations
Globals1D;

%% Order of polymomials used for approximation 
N = 4;

%% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,10);

%% Initialize solver and construct grid and metric
StartUp1D;

%% Set initial conditions
u0 = 0.1 + 0.5*exp(-20*x.^2);

%% Solve Problem
FinalTime = 1.0;
[u] = Advec1D(u0,FinalTime);

%% plot
%plot(x,u);