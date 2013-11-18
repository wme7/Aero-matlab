%% SQUARESTOKE Stokes equations on the unit square
%
%   SQUARESTOKE computes P2-P1 approximations of the Stokes equations in
%   the unit square on a sequence of meshes obtained by uniform refinement.
%   It plots the approximation error (pressue in L2 norm and velocity in H1
%   norm) vs the number of dof.
% 
% See also StokesP2P1, collidingflow, squarePoisson 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Set up
maxIt = 4;
N = zeros(maxIt,1); erru = zeros(maxIt,1); errp = zeros(maxIt,1);

%% Generate initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet');

%% PDE and options
pde = Stokesdata2;

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [u,p,edge,A] = StokesP2P1(node,elem,pde,bdFlag);
    N(k) = length(u)+length(p);
    if N(k) < 2e3 % show solution for small size
        figure(1);  showresult(node,elem,p);    
    end
    % compute error
    uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    erru(k) = sqrt((u-uI(:))'*A*(u-uI(:)));
    errp(k) = getL2error(node,elem,pde.exactp,p);
    % refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end

%% Plot convergence rates
figure(2); clf;
showrate2(N,erru,2,'-*','||Du_I-Du_h||',N,errp,2,'k-+','|| p - p_h||');