%% Exampels of Stokes Equations
%
% We test several finite element discretizations for Stokes equations.

%% Different elements
% * P2-P1 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'P2P1';
StokesFEM;
%%
% * isoP2-P1 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'isoP2P1';
StokesFEM;
%%
% * P2-P0 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'P2P0';
StokesFEM;
%%
% * isoP2-P0 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'isoP2P0';
StokesFEM;
%%
% * CR element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'CR';
StokesFEM;
%%
% * CRP1 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'CRP1';
StokesFEM;
%%
% * Mini element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'Mini';
StokesFEM;
%%
% * P1bP1 element
clear all
pde = Stokesdata1;
[node,elem] = squaremesh([-1 1 -1 1], 0.125);
option.fem = 'P1bP1';
StokesFEM;

%% Different solvers
% option.solver = 

% %% Different data and boundary conditions
% % Poiseuille flow in a square domain
% clear all
% pde = Stokesdata3;
% [node,elem] = squaremesh([-1 1 -1 1], 1);
% bdFlag = setboundary(node,elem,'Dirichlet','~(x==1)','Neumann','x==1');
% % option.fem = 'P2P1';
% % option.fem = 'isoP2P1';
% % option.fem = 'P2P0';
% % option.fem = 'isoP2P0';
% % option.fem = 'CR';
% % option.fem = 'Mini';
% option.fem = 'P1bP1';
% StokesFEM;
% 
% %% Different mesh