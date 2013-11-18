%% RATE OF CONVERGENCE OF LINEAR FINITE ELEMENT METHOD
%
% This example is to show the rate of convergence of linear finite element
% approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the following boundary condition:
%
% # Non-empty Dirichlet boundary condition. $u=g_D \hbox{ on }\Gamma_D, \quad \nabla u\cdot n=g_N \hbox{ on }\Gamma_N. \Gamma _D = \{(x,y): x=0, y\in [0,1]\}, \; \Gamma _N = \partial \Omega \backslash \Gamma _D$. 
% # Pure Neumann boundary condition. $\Gamma _N = \partial \Omega$.
% # Robin boundary condition. $g_R u + \nabla u\cdot n=g_N \hbox{ on }\partial \Omega$

%% 
clear all; close all; clc;
[node,elem] = squaremesh([0,1,0,1],0.25); 
option.L0 = 2;
option.maxIt = 4;
option.maxN = 1e6;
option.printlevel = 1;
option.plotflag = 0;
option.elemType = 'P3';

%% Non-empty Dirichlet boundary condition.
pde = sincosdata;
bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
% bdFlag = setboundary(node,elem,'Dirichlet');
femPoisson(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
pde = sincosNeumanndata;
bdFlag = setboundary(node,elem,'Neumann');
femPoisson(node,elem,pde,bdFlag,option);

%% Pure Robin boundary condition.
option.plotflag = 0;
pdeRobin = sincosRobindata;
bdFlag = setboundary(node,elem,'Robin');
femPoisson(node,elem,pdeRobin,bdFlag,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (2nd order) and L2-norm
% (3rd order) is observed. The order of ||DuI-Duh|| is almost 3rd order and
% thus superconvergence exists.
%
% MGCG converges uniformly in all cases.
