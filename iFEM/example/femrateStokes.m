%% RATE OF CONVERGENCE OF FINITE ELEMENT METHOD FOR THE STOKES EQNS IN 2D
%
% This example is to show the rate of convergence of various finite element
% approximation of the Stokes equation on the unit square:
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%                              u = g_D  on \Gamma.
%
% with the pure Dirichlet boundary condition.
%
% Written by Ming Wang. Modified by Long Chen.

clear all; close all;
%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
% pde = Stokesdata1;
pde = Stokesdata2;
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 2;
option.maxIt = 3;
option.printlevel = 1;
option.plotflag = 1;

%% P2-P1 element
option.elemType = 'P2P0';
option.solver = 'asmg';
femStokes(node,elem,pde,bdFlag,option);

%% P2-P1 element
option.elemType = 'P2P1';
option.solver = 'mg';
femStokes(node,elem,pde,bdFlag,option);

% %% P2-P0 element
% option.fem = 'P2P0';
% option.solver = 'mg';
% femStokes(node,elem,pde,bdFlag,option);

% option.fem = 'P2P0';
% option.fem = 'isoP2P1';
% option.fem = 'isoP2P0';
% option.fem = 'CR';
% option.fem = 'CRP1';  
% option.fem = 'Mini';
% option.fem = 'P1bP1';

%% Conclusion
% For the Stokesdata3, the exact solution is polynoimal, thus FEM solution 
% is exact.
%  
% The solution of CRP1 for mixed boundary condition problem is not right.
% 