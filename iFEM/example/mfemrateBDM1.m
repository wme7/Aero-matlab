%% RATE OF CONVERGENCE OF MIXED FINITE ELEMENT METHOD
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
[node,elem] = squaremesh([0,1,0,1],0.25); 
pde = mixBCdata;
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'BDM1-P0';
option.solver = 'uzawapcg';
% option.solver = 'dmg';
% option.solver = 'direct';

%% Non-empty Dirichlet boundary condition.
option.plotflag = 1;
bdFlag = setboundary(node,elem,'Dirichlet');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
option.plotflag = 0;
bdFlag = setboundary(node,elem,'Neumann');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Mix Dirichlet and Neumann boundary condition.
option.plotflag = 0;
% option.solver = 'uzawapcg';
bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Conclusion
%
% The optimal rates of convergence for u and sigma are observed, namely,
% 1st order for L2 norm of u, H(div) norm of sigma, and 2nd order for 
% L2 norm of sigma. 
%
% The 2nd order convergent rates between two discrete functions ||uI-uh|| 
% is known as superconvergence.
%
% The 3rd order convergent rates between two discrete functions 
% ||sigmaI-sigmah|| for the Pure Neumman boundary condition 
% is known as superconvergence. The 2nd accurate numerical quadrature
% is required for the integral of rhs to observe such superconvergence.
%
% dmg and uzawapcg converges uniformly in all cases. 
% The tolerence should be small enough for Distributive MG to achieve the 
% accurate solution. 