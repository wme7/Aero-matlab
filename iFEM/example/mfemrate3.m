%% RATE OF CONVERGENCE OF MIXED FINITE ELEMENT METHOD IN 3D
%
% This example is to show the rate of convergence of mixed finite element
% (RT0-P0) approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the following boundary condition:
%
% # Pure Dirichlet boundary condition. $\Gamma _D = \partial \Omega$. 
% # Pure Neumann boundary condition. $\Gamma _N = \partial \Omega$.
% # Mix Dirichlet and Neumann boundary condition. $u=g_D \hbox{ on }\Gamma_D, 
% \quad \nabla u\cdot n=g_N \hbox{ on }\Gamma_N. \Gamma _N = \{(x,y): x=0, 
% y\in [0,1]\}, \; \Gamma _D = \partial \Omega \backslash \Gamma _N$. 
%
% Written by Ming Wang.

%% 
clear all; close all;
% [node,elem] = cubemesh([-1,1,-1,1,-1,1],2);
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];  % nodes
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7]; % elements
pde = mixBCdata3;
option.L0 = 2;
option.maxIt = 3;
option.printlevel = 1;
option.elemType = 'RT0-P0';
option.refType = 'bisect';
% option.solver = 'uzawapcg';
% option.solver = 'dmg';
option.solver = 'direct';

%% Pure Dirichlet boundary condition.
% option.plotflag = 0;
% bdFlag = setboundary3(node,elem,'Dirichlet);
bdFlag = setboundary3(node,elem,'Neumann');
mfemPoisson3(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
% option.plotflag = 0;
bdFlag = setboundary(node,elem,'Neumann');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Mix Dirichlet and Neumann boundary condition.
option.plotflag = 1;
option.solver = 'uzawapcg';
bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Conclusion
%
% The optimal rates of convergence for u and sigma are observed, namely,
% 1st order for L2 norm of u, L2 norm of sigma and H(div) norm of sigma. 
% The 2nd order convergent rates between two discrete functions ||uI-uh|| 
% and ||sigmaI-sigmah|| are known as superconvergence.
%
% dmg and uzawapcg converges uniformly in all cases. Distributive MG is two
% times faster than uzawapcg.