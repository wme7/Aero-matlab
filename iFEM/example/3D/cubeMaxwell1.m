%% CUBEMAXWELL11 solves Maxwell type equations in a cube using linear element.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;

%% Defacult setting
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
pde = Maxwelldata2;
% pde = planewavedata1;
% bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary3(node,elem,'Dirichlet');
option.solver = 'cg';

%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
%     [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
    % solve the equation
    [u,edge,eqn] = Maxwell1(node,elem,HB,pde,bdFlag,option); 
    % compute error
    energyErr(k) = getHcurlerror3NE1(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3NE1(node,elem,pde.exactu,real(u));
    N(k) = length(u);
end

%% Plot convergence rates
N = N(1:k); energyErr = energyErr(1:k); L2Err = L2Err(1:k);
figure(1); clf; 
showrate2(N,energyErr,1,'r-+','||u-u_h||_A',...
          N,L2Err,1,'b-+','||u-u_h||');