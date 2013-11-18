%% MULTIGRID OF FOR THE STOKES EQNS IN 2D
%
% This example is to show the convergence of multigrid methods for various
% finite element approximation of the Stokes equation on the unit square:
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%                              u = g_D  on \Gamma.
%
% with the pure Dirichlet boundary condition.
%

%% Setting
clear all; close all;
[node,elem] = squaremesh([0,1,0,1],0.5);
% [node,elem] = circlemesh(0,0,1,0.25);
[node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
pde = Stokesdata1; 
%  pde = StokesZulehnerdata;
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 1;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 0;
option.rateflag = 0;

%% MG options
option.solver = 'mg';
option.smoothingstep = 2;
option.smootherbarSp = 'SGS';

%% RT0-P0
display('RT0-P0')
option.elemType = 'RT0-P0';
option.refType = 'bisect';
femStokesHdiv(node,elem,pde,bdFlag,option);

% %% RT0-P0
% display('BDM1B-P0')
% option.elemType = 'BDM1B-P0';
% femStokesHdiv(node,elem,pde,bdFlag,option);

%% CR-P0 element
display('CR-P0')
option.elemType = 'CRP0';
femStokes(node,elem,pde,bdFlag,option);

%% P2-P0 element
display('P2-P0')
option.elemType = 'P2P0';
femStokes(node,elem,pde,bdFlag,option);

%% isoP2-P0 element
display('isoP2-P0')
option.elemType = 'isoP2P0';
femStokes(node,elem,pde,bdFlag,option);

%% isoP2-P1 element
display('isoP2-P1')
option.elemType = 'isoP2P1';
femStokes(node,elem,pde,bdFlag,option);

% %% P1b-P1 element
% display('P1b-P0')
% option.elemType = 'P1bP1';
% femStokes(node,elem,pde,bdFlag,option);

%% P2-P1 element
display('P2-P1')
option.elemType = 'P2P1';
% option.smoothingStep = 3;
% option.smootherbarSp   = 'VCYCLE';
% option.smootherbarSpPara = 0.75;
femStokes(node,elem,pde,bdFlag,option);
