%% Driver script for solving the 1D advection equations
Globals1D;

%% Order of polymomials used for approximation 
N = 5;

%% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,4.0,10);

%% Initialize solver and construct grid and metric
StartUp1D;

%% Set initial conditions
u = sin(x);

% u_0 = zeros(1,K*Np);
% jump = [0.9 1.1];
% for i = 1:K*Np
%     if jump(1) <= x(i) && x(i) <= jump(2)
%         u_0(i) = 2;
%     else
%         u_0(i) = 0;
%     end
% end
% Load Initial Condition
u = reshape(u_0,Np,K);

%% Solve Problem
FinalTime = 1.2;
[u] = Advec1D(u,FinalTime);
