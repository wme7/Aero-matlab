% Driver script for solving the 1D Burgers equations
Globals1D;

% Orer of polymomials used for approximation 
N = 6;

% Generate simple mesh
xL = -1.0; xR = 1.0;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
epsilon = 0.1;
%u = -tanh((x+0.5)/(2*epsilon)) + 1.0;
%u=2*sin(x);
xmid=(x(end)+x(1))/2; x1=x<xmid; x2=x>=xmid; u=1*x1+0.1*x2;

% Solve Problem
FinalTime = 0.15;
[u] = Burgers1D(u,epsilon,xL,xR,FinalTime);