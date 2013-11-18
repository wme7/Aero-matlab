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
uIuhErr = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
%     [node,elem,bdFlag,HB] = uniformbisect3(node,elem,bdFlag,HB);
    % solve the equation
    if(abs(pde.epsilon)>1.0e-8)
    [u,edge,eqn] = Maxwell(node,elem,HB,pde,bdFlag,option); 
    else
    [u,edge,eqn] = Maxwellsaddle(node,elem,HB,pde,bdFlag,option); 
    end
    % compute error
    tic;
    energyErr(k) = getHcurlerror3NE(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3NE(node,elem,pde.exactu,real(u));
    uI = edgeinterpolate(pde.exactu,node,edge);
    uIuhErr(k) = sqrt((real(u)-uI)'*(eqn.A)*(real(u)-uI));        
    time = toc;
    fprintf('Time to compute the error %4.2g s\n',time);
    N(k) = length(u);
end

%% Plot convergence rates
figure(1); clf; 
showrate3(N,energyErr,1,'r-+','||u-u_h||_A',...
          N,L2Err,1,'b-+','||u-u_h||',...
          N,uIuhErr,1,'m-+','||uI-u_h||_A');