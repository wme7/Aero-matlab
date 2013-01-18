% Driver script for solving the 1D Euler equations
Globals1D;

% Polynomial order used for approximation 
N = 6;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0, 1.0, 100);

% Initialize solver and construct grid and metric
StartUp1D;

% Set up initial conditions -- Sod's problem
MassMatrix = inv(V')/V;
cx = ones(Np,1)*sum(MassMatrix*x,1)/2; 

%u = ones(Np,K).*( sin(2*pi*cx) );
u = sin(2*pi*x);
FinalTime = 5.0;

% Solve Problem
[u] = iBurgers1D(u,FinalTime);
snapnow;
