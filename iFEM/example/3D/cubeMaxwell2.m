%% CUBEMAXWELL solves Maxwell type equations in a cube using lowest order element.
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
%     [node,elem,bdFlag,HB] = bisect3(node,elem,'all',bdFlag,HB);
    [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
    % solve the equation
    [u,T,eqn] = Maxwell2(node,elem,HB,pde,bdFlag,option); 
    % compute error
    energyErr(k) = getHcurlerror3NE2(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3NE2(node,elem,pde.exactu,real(u));
%     uI = edgeinterpolate2(pde.exactu,node,T.edge,T.face,T.face2edge);
%     L2Err(k) = sqrt(abs((u-uI)'*eqn.M*(u-uI)));
%     energyErr(k) = sqrt((u-uI)'*eqn.A*(u-uI) + L2Err(k)^2);
    N(k) = length(u);
end

%% Plot convergence rates
N = N(1:k); energyErr = energyErr(1:k); L2Err = L2Err(1:k);
figure(1); clf; 
showrate2(N,energyErr,1,'r-+','||u-u_h||_A',...
          N,L2Err,1,'b-+','||u-u_h||');