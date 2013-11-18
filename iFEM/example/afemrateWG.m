%% RATE OF CONVERGENCE OF ADAPTIVE FINITE ELEMENT METHOD USING WG ELEMENT
%
% This example is to show the rate of convergence of the lowest order
% finite element approximation of the second order elliptic equation.
%
% # Lshape problem.
% # Kellogg problem.

clear all; close all;
%% Kellogg problem
[node,elem] = squaremesh([-1 1 -1 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Kelloggdata;
% singular node
distance = (node(:,1)-0).^2 + (node(:,2)-0).^2;
[mindist,singularnode] = min(distance); %#ok<ASGLU>
pde.singularnode = singularnode;
option.L0 = 0;
option.maxIt = 500;
option.maxN = 1e4;
option.theta = 0.2;
option.elemType = 'WG';
option.plotflag = 1;
option.printlevel = 0;
[err,time,solver,eqn,node,elem] = afemPoisson(node,elem,pde,bdFlag,option);
figure;
showrate2(err.N,err.H1,270,'k-.','||Du-Du_h||',err.N,err.eta,250,'k-','eta');
latexerrtable(err.N,[err.H1 err.eta])
% rate = -log(err.H1)./log(err.N);
% figure; plot(rate,'r-.');
% 
% %% Lshape problem
% [node,elem] = squaremesh([-1,1,-1,1],1);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');
% bdFlag = setboundary(node,elem,'Dirichlet');
% pde = Lshapedata;
% format shorte
% option.L0 = 1;
% option.maxIt = 50;
% option.printlevel = 1;
% option.elemType = 'WG';
% option.plotflag = 1;
% err = afemPoisson(node,elem,pde,bdFlag,option);

