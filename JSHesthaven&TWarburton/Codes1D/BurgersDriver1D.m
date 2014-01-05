% Driver script for solving the 1D Burgers equations
Globals1D;

% Order of polymomials used for approximation 
N = 8;

% Generate simple mesh
xL = -0.0; xR = 2.0*pi;
[Nv, VX, K, EToV] = MeshGen1D(xL,xR,20);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
epsilon = 0.1;
%u = -tanh((x+0.5)/(2*epsilon)) + 1.0;
u=2*sin(x);
%xmid=(x(end)-x(1))/2; x1=x<xmid; x2=x>=xmid; u=1*x1+0*x2;

% Solve Problem
FinalTime = 2.0;%0.75;
[u] = Burgers1D(u,epsilon,xL,xR,FinalTime);

plot(x,u,'-+'); grid on;