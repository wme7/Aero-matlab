%% AMG TEST IV: Three Dimensional Problem
%
% We consider the linear finite element discretization of Poisson equation
% in three dimensions with Dirichlet boundary conditions and compare
% geometric multigrid and algebraic multigrid.

%% Conclusion
%
% The algebraic multigrid is robust to the size of the matrix. The
% iteration steps increases slightly. The preformance is better on
% structure grids than unstructure grids. The sparsity of the coarse grid
% is increased.

%% Uniform Mesh: Geometric Multigrid
clear all
option.solver = 'mg';
cubePoisson;

%% Uniform Mesh: Algebraic Multigrid
clear all
option.solver = 'amg';
cubePoisson;

%% Unstructured Mesh
clear all;
load oilpump;
showboundary3(node,elem);
option.solver = 'amg';
cubePoisson;