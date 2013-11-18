%% Edge Element Discretization of Maxwell Equations: Linear Element
% We test Maxwell solvers in iFEM.

clear all; close all
rowNames ={'h=1/2';'h=1/4';'h=1/8'};
colHeaders = {'H^curl Error','L^2 Error'};

%% The data of the pde
%
% * pde = Maxwelldata1; % zero Neumann boundary condition and curl u = 0
% * pde = Maxwelldata2; % non-homogenous Neumann boundary condition
% * pde = Maxwelldata3; % polynomial data and curl u = 0
% * pde = Maxwelldata4; % zero Dirichlet boundary condition
% * pde = Maxwelldata5; % linear polynomial data
% * pde = planewavedataC; % plane wave with complex coefficients
% * pde = planewavedata1; % plane wave with real coefficients


%% Positive Definite Case
% curl curl E + E = f.
option.solver = 'cg';

help Maxwelldata2
%% 
% Dirichlet boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = Maxwelldata2;
bdFlag = setboundary3(node,elem,'Dirichlet');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%% 
% Neumann boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = Maxwelldata2;
bdFlag = setboundary3(node,elem,'Neumann');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%%
% Optimal first order of convergence is achieved. HX preconditioned CG
% converges around 20 steps. 

%% Indefinite with real coefficients
% curl curl E - E = f.

help planewavedata1
%% 
% Dirichlet boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = planewavedata1;
bdFlag = setboundary3(node,elem,'Dirichlet');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%% 
% Neumann boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = planewavedata1;
bdFlag = setboundary3(node,elem,'Neumann');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%%
% Optimal first order of convergence is achieved. HX preconditioned CG
% converges around 40 steps although the system is indefinite.

%% Indefinite: complex coefficients, real solution
% curl curl E - (1-i)E = f.
option.solver = 'gmres';

help planewavedataC
%% 
% Dirichlet boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = planewavedataC;
bdFlag = setboundary3(node,elem,'Dirichlet');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);


%% 
% Neumann boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = planewavedataC;
bdFlag = setboundary3(node,elem,'Neumann');
cubeMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%%
% Optimal first order of convergence in energy norm is achieved. But L2
% norm is not optimal. HX preconditioned GMRES converges around 90 steps
% although the system is indefinite.
%
% For Neumann boundary condition, the rate of L2 error is only half order
% but the solver converges in less steps (40 - 50 steps).

%% Indefinite: real coefficents, complex solution
option.solver = 'cg';

help planewavedata
%%
% Dirichlet boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');
planewaveMaxwell1;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%% 
% Neumann boundary condition
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
bdFlag = setboundary3(node,elem,'Neumann');
planewaveMaxwell;
%%
makeHtmlTable([energyErr L2Err],[],rowNames,colHeaders);

%%
% The computation of error can't handle complex functions. So in
% planewaveMaxwell the error between uI and uh is computed. Therefore
% slightly better rate of convergence is observed. The rate of L2 error is
% almost second order.