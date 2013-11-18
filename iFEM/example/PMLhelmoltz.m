%% RATE OF CONVERGENCE OF LINEAR FINITE ELEMENT METHOD FOR HELMHOLTZ EQUATION
%
% This example is to show the rate of convergence of linear finite element
% approximation of the PML Helmholtz equation on the unit square:
%
% for the following boundary condition:
%

%% 
clear all; close all;
global k        % wave number
h = 1/100;
[node] = squaremesh([0,1,0,1],h);   %intresting domain
[node,elem]  =  PMLdomain(node,h);  %PML domain
pde = Helmholtzdata3;
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;

%% Dirichlet boundary condition.
option.plotflag = 1;
bdFlag = setboundary(node,elem,'Dirichlet');
option.solver = 'direct';
%for k = 1:2:3
k=7; 
u = HelmholtzPML(node,elem,pde,bdFlag,option);
showsolution(node,elem,u);
    %%
%end


%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
%
% MGCG converges uniformly in all cases.