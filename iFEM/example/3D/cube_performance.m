%% CUBE_PERFORMANCE solves Poisson equation in a cube.
% 
% It tests the performance of iFEM for large number of unknowns.
%
% See also  cubeAFEM, Lshape
%
% <a href="matlab:ifemdoc cube_performance">iFEMdoc cube_performance</a>
%
% Copyright (C) Long Chen.

close all; clear all;
%% Parameters
maxIt = 4; N = zeros(maxIt,1); 

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);

%% Set up boundary condition
bdFlag = setboundary3(node,elem,'Dirichlet');

%%  Get a fine mesh by uniform bisection
for k = 1:2
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
end

%% Get the data of the pde
pde = sincosdata3;

%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
N = zeros(maxIt,1);     cost = zeros(maxIt,1);
for k = 1:maxIt
    tic;
    u = Poisson3(node,elem,pde,HB,bdFlag);       % Solve
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);   % Refine
    cost(k) = toc;
    N(k) = length(u); 
end

%% Plot computational cost
figure; 
showrate(N,cost,4,'-*','cpu time');
title('Computational cost vs Number of dof', 'FontSize', 14);
display(['DOF     ' '     CPU (second)']); 
display(num2str([N cost]));