%% CUBEMAXWELL solves Maxwell type equations in a cube.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Parameters
maxIt = 4; 
N = zeros(maxIt,1); 
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);

%% Set up boundary condition
% bdFlag = setboundary3(node,elem,'Dirichlet');
bdFlag = setboundary3(node,elem,'Neumann');
% showmesh3(node,elem); view([130,28]);

%%  Get a fine mesh by uniform bisection
for k = 1:1
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
%     [node,elem,HB,bdFlag] = bisect3(node,elem,'all',HB,bdFlag);
end

%% Get the data of the pde
 %pde = Maxwelldata1; % zero Neumann boundary condition and curl u = 0
% pde = Maxwelldata2; % non-homogenous Neumann boundary condition
% pde = Maxwelldata3; % polynomial data and curl u = 0
% pde = Maxwelldata4; % zero Dirichlet boundary condition
%pde = Maxwelldata5; % linear polynomial data
pde = planewavedataH; % plane wave

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
%     [node,elem,HB,bdFlag] = bisect3(node,elem,'all',HB,bdFlag);
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
    % solve the equation
%     [u,edge,A,M] = Maxwell(node,elem,HB,pde,bdFlag);
%     option.solver = 'minres'; 
    option.printlevel = 2;
    [u,edge,eqn] = MaxwellH(node,elem,HB,pde,bdFlag,option); 
    % compute uI by two-points quadrature
%     uI = edgeinterpolate(pde.exactu,node,edge);
%     energyErr(k) = sqrt((u-uI)'*(eqn.A-eqn.M)*(u-uI));
%     L2Err(k) = sqrt(-(u-uI)'*eqn.M*(u-uI));
    energyErr(k) = getHcurlerror3NE(node,elem,pde.curlu,u);
    L2Err(k) = getL2error3NE(node,elem,pde.exactu,u);
    N(k) = length(u);
end

%% Plot convergence rates
N = N(1:k); 
figure(1); clf; 
r1 = showrate(N(1:k),energyErr(1:k),2,'r-+');
hold on
r2 = showrate(N(1:k),L2Err(1:k),2,'b-+');
legend('||u-u_h||_A',['N^{' num2str(r1) '}'],...
       '||u-u_h||',['N^{' num2str(r2) '}'],'LOCATION','Best');
