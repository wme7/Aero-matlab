% Driver script for solving the 1D Euler equations
Globals1D;

% Polynomial order used for approximation 
N = 4;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0, 2*pi, 80);

% Initialize solver and construct grid and metric
StartUp1D;

% Set up initial conditions -- Sod's problem
MassMatrix = inv(V')/V;
%cx = ones(Np,1)*sum(MassMatrix*x,1)/2; 

%u = ones(Np,K).*( sin(2*pi*cx) );
u0 = 0.5 + sin(x);
FinalTime = 1.5;

% Solve Problem
[u] = iBurgers1D(u0,FinalTime);
snapnow;

% Plot Figure
plotrange=[0,2*pi,min(min(u0))-0.2,max(max(u0))+0.2];
plot(x,u,x,u0,'-+'); axis(plotrange); ylabel('u(x)'); xlabel('x');