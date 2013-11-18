%% COMPARISON OF CONVERGENCE RATE OF VARIOUS FINITE ELEMENT METHODS
%
% This example is to compare the rate of convergence of various finite
% element approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the Dirichlet boundary condition. Results holds for other boundary
% conditions. 

%% Problem
clear all; close all; clc;
[node,elem] = squaremesh([0,1,0,1],0.25); 
pde = sincosdata;
bdFlag = setboundary(node,elem,'Dirichlet');
% pde = sincosNeumanndata;
% bdFlag = setboundary(node,elem,'Neumann');
option.L0 = 2;
option.maxIt = 5;
option.printlevel = 1;
option.plotflag = 0;

%% P1 element 
err = femPoisson(node,elem,pde,bdFlag,option);

%%
% The optimal rate of convergence of the H1-norm (1st order), L2-norm (2nd
% order), and discrte maximum norm is observed. The 2nd order convergent
% rate between two discrete functions ||DuI-Duh|| is known as
% superconvergence.

%% CR nonconforming P1 element
option.elemType = 'CR';
err = femPoisson(node,elem,pde,bdFlag,option);

%%
% The optimal rate of convergence of the H1-norm (1st order), L2-norm (2nd
% order), and discrte maximum norm is observed. But no superconvergence of
% ||DuI-Duh||.

%% P2 element
option.elemType = 'P2';
err = femPoisson(node,elem,pde,bdFlag,option);

%% P3 element
option.elemType = 'P3';
err = femPoisson(node,elem,pde,bdFlag,option);

%% Conclusion

