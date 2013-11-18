%% RATE OF CONVERGENCE OF WEAK GALERKIN FINITE ELEMENT METHOD IN 3-D
%
% This example is to show the rate of convergence of weak Galerkin finite
% element approximation of the Poisson equation on the unit cube:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^3$$
%
% for the following boundary condition:
%
% # Non-empty Dirichlet boundary condition. $u=g_D \hbox{ on }\Gamma_D, \quad \nabla u\cdot n=g_N \hbox{ on }\Gamma_N. \Gamma _D = \{(x,y): x=0, y\in [0,1]\}, \; \Gamma _N = \partial \Omega \backslash \Gamma _D$. 
% # Pure Neumann boundary condition. $\Gamma _N = \partial \Omega$.
% # Robin boundary condition. $g_R u + \nabla u\cdot n=g_N \hbox{ on }\partial \Omega$

%% 
clear all; close all;
[node,elem] = cubemesh([0,1,0,1,0,1],0.25); 
pde = sincosdata3;
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'WG';
option.plotflag = 0;

%% Non-empty Dirichlet boundary condition.
bdFlag = setboundary3(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
% bdFlag = setboundary3(node,elem,'Dirichlet');
femPoisson3(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
% pde = sincosNeumanndata;
bdFlag = setboundary3(node,elem,'Neumann');
femPoisson3(node,elem,pde,bdFlag,option);

%% Pure Robin boundary condition.
pdeRobin = sincosRobindata3;
bdFlag = setboundary3(node,elem,'Robin');
femPoisson3(node,elem,pdeRobin,bdFlag,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. No superconvergence for ||DuI-Duh||.
%
% MGCG converges uniformly in all cases.