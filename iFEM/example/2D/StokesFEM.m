%% STOKEFEM Finite Element Methods for Stokes Equations
%
% 
% See also StokesP2P1, collidingflow, squarePoisson 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

clear all
close all
%% Set up
if ~exist('maxIt','var'), maxIt = 4; end
Ndof = zeros(maxIt,1); erru = zeros(maxIt,1); errp = zeros(maxIt,1);

%% Generate initial mesh
if ~exist('node','var')
    [node,elem] = squaremesh([0 1 0 1], 0.125);
end

%% Boundary conditions
if ~exist('bdFlag','var')
    bdFlag = setboundary(node,elem,'Dirichlet');
end

%% PDE
if ~exist('pde','var')
    pde = Stokesdata2;
end

%% Elements
if ~exist('option','var')
    option.elemType = 'P2P1';
end

%% Solver
if ~isfield(option,'solver')
    option.solver = 'direct';
    if (strcmp(option.elemType,'P2P1') || strcmp(option.elemType,'P2P0')|| strcmp(option.elemType,'CRP0')|| strcmp(option.elemType,'CRP1'));
        option.solver = 'mg';
    end
end
option = femStokesoption(option);
%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [u,p,edge,A] = Stokes(node,elem,pde,bdFlag,option);

    Ndof(k) = length(u)+length(p);
%     if N(k) < 2e3 % show solution for small size
%         figure(1);  showresult(node,elem,p);    
%     end
    % compute error
    if strcmp(option.elemType,'P2P0') || strcmp(option.elemType,'P2P1') || ...
       strcmp(option.elemType,'isoP2P0') || strcmp(option.elemType,'isoP2P1')
        uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
    elseif strcmp(option.elemType,'CRP0') || strcmp(option.elemType,'CRP1')
        uI = pde.exactu((node(edge(:,1),:)+node(edge(:,2),:))/2);
    elseif strcmp(option.elemType,'MINI') 
        uI = pde.exactu(node);
    elseif strcmp(option.elemType,'P1bP1')
        % bubble part won't be taken in the error computation
        N = size(node,1); NT = size(elem,1);
        u0 = pde.exactu(node);
        uI = u;
        uI([(1:N)'; NT+N+(1:N)']) = u0(:);
    end
    erru(k) = sqrt((u-uI(:))'*A*(u-uI(:)));
    errp(k) = getL2error(node,elem,pde.exactp,p);
    % refine mesh
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end

%% Plot convergence rates
figure(2); clf;
showrate2(Ndof,erru,2,'-*','||Du_I-Du_h||',Ndof,errp,2,'k-+','|| p - p_h||');