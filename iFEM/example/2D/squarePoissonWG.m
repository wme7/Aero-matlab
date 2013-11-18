%% SQUAREPOISSON Poisson equation in a square domain.
%
%   squarePoisson computes approximations of the Poisson equation in the
%   unit square on a sequence of meshes obtained by uniform refinement. It
%   plots the approximation error (in L2 norm or H1 norm) vs the number of
%   nodes.
% 
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Parameters 
maxIt = 4; N = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = squaremesh([0 1 0 1], 0.25);
%bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%bdFlag = setboundary(node,elem,'Dirichlet');
for k = 1:3
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Get the data of the pde
pde = sincosdata;
% pde = mixBCdata;
option.solver = 'mg';

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [u,Du,eqn] = PoissonWG(node,elem,pde,bdFlag,option);
    N(k) = size(elem,1);
    % compute error
    uI = zeros(N(k)+size(eqn.edge,1),1);
    uI(1:N(k)) = pde.exactu((node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3);
    uI(N(k)+1:end) = pde.exactu((node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2);
    erruIuh(k) = sqrt((u-uI)'*eqn.A*(u-uI));
    % refine mesh
   [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Plot convergence rates
figure(2);
showrate(N,erruIuh);
%% Error table
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
