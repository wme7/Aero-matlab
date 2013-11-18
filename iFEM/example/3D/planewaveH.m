%% PLANEWAVE plane wave solutions to Maxwell equations in a cube.
%
% solve for H
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
energyErrImag = zeros(maxIt,1);
L2ErrImag = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);

%% Set up boundary condition
% bdFlag = setboundary3(node,elem,'Dirichlet');
bdFlag = setboundary3(node,elem,'Neumann');
% showmesh3(node,elem); view([130,28]);
% findnode3(node);

%% Get the data of the pde
global d P omega
d = [1, pi/2, pi/2]; % sphereical coordinate
r = d(1); theta = d(2); phi = d(3);
d = [r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)]; % Cartesain coordinate
P = [1, 0, 0];
omega = 1;
pde = planewavedataH; % plane wave
pde.omega = omega;

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
%     [node,elem,HB,bdFlag] = bisect3(node,elem,'all',HB,bdFlag);
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
    % solve the equation
    option.printlevel = 1;
    [u,edge,eqn] = MaxwellH(node,elem,HB,pde,bdFlag,option);
    uI = edgeinterpolate(pde.exactu,node,edge);
%     [u,edge,eqn] = Maxwell1H(node,elem,HB,pde,bdFlag,option);
%     uI = edgeinterpolate1(pde.exactu,node,edge);
%     [u,T,eqn] = Maxwell2H(node,elem,HB,pde,bdFlag); 
%     uI = edgeinterpolate2(pde.exactu,node,T.edge,T.face,T.face2edge);
    % compute uI by two-points quadrature
    L2Err(k) = sqrt(abs(real(u-uI)'*eqn.M*real(u-uI)));    
    energyErr(k) = sqrt(abs(real(u-uI)'*eqn.A*real(u-uI)) + L2Err(k)^2);
    L2ErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.M*imag(u-uI)));
    energyErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.A*imag(u-uI)) + L2ErrImag(k)^2);
%     energyErr(k) = getHcurlerror3NE(node,elem,pde.curlu,u);
%     L2Err(k) = getL2error3NE(node,elem,pde.exactu,u);
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
figure(2); clf; 
r3 = showrate(N(1:k),energyErrImag(1:k),2,'r-+');
hold on
r4 = showrate(N(1:k),L2ErrImag(1:k),2,'b-+');
legend('||u-u_h||_A',['N^{' num2str(r3) '}'],...
       '||u-u_h||',['N^{' num2str(r4) '}'],'LOCATION','Best');
