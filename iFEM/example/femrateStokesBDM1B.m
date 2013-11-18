%% RATE OF CONVERGENCE OF THE TMAC DISCRETIZATION FOR THE STOKES EQUANTIONS
%
% This example is to show the rate of convergence of (BDM1+bubble)-P0 
% finite element approximation of the Stokes equation on the unit square:
%
% -------------------------------------------------------------------------
%
% - grad div u + curl rot u + grad p  = f   in $\Omega$      
%
% - div u   = 0   in $\Omega$
%
% -------------------------------------------------------------------------
%
% for the Dirichlet boundary condition:
%
% -------------------------------------------------------------------------
%
% u $\cdot$ t   = $g_t$   on $\Gamma$   
%
% u $\cdot$ n   = $g_n$   on $\Gamma$   
%
% -------------------------------------------------------------------------
%
% Test example: [0,1] $\times$ [0,1]
%
%   u = -(2^8)*(x.^2-2*x.^3+x.^4).*(2*y-6*y.^2+4*y.^3)
%   v = 2^8*(2*x-6*x.^2+4*x.^3).*(y.^2-2*y.^3+y.^4)
%   p = -(2^8)*(2-12*x+12*x.^2).*(y.^2-2*y.^3+y.^4)
% 
% Created by Ming Wang, at Nov., 2011.

%%
clear all; close all;
[node,elem] = squaremesh([0,1,0,1],0.5);
% load('BraessData.mat');
% pde = testdata;
%  pde = StokesNHsincos;
%  pde = StokesBEwup;
% pde = StokesZulehnerdata;
%  pde = StokesLChenProjectdata;
pde = Stokesdata1;
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 2;
option.elemType = 'BDM1B-P0';
option.maxIt = 4;
option.solver = 'mg';
option.smoothingstep = 2;
option.smootherbarSp = 'SGS';
% option.printlevel = 3;
% %% Refine Strategy: uniform bisect
% option.refType = 'red';
option.refType = 'bisect';
femStokesHdiv(node,elem,pde,bdFlag,option);

% %% Refine Strategy: uniform refine
% close all;
% option.L0 = 3;
% % [node,elem]=squaremesh([0,1,0,1],1);
% % bdFlag = setboundary(node,elem,'Dirichlet');
% option.plotflag = 1;
% option.refType = 'red';
% femStokesHdiv(node,elem,pde,bdFlag,option);

%% Conclusion
% The optimal rate of convergence of L2 norm (1st order) is observed for u.
% The optimal rate of convergence of L2 norm (1st order) is observed for p.
%
% When we use refine the mesh uniformly, 2st superconverge rate is observed
% for the differecne of u_h and u_I in L2 norm.