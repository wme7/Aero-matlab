%% PLANEWAVEMAXWELL plane wave solutions to Maxwell equations in a cube.
%
% solve for E
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
energyErrImag = zeros(maxIt,1);
L2ErrImag = zeros(maxIt,1);

%% Get the data of the pde
global d P omega
d = [1, pi/2, pi/2];  % in sphereical coordinate
P = [1, 0, 0];
omega = 1;
r = d(1); theta = d(2); phi = d(3);
d = [r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)];
% pde = planewavedataC; % plane wave with complex coefficients
pde = planewavedata; % plane wave with real coefficients and complex solution

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);
    % solve the equation
%     [u,edge,A,M] = Maxwell(node,elem,HB,pde,bdFlag);
    [u,edge,eqn] = Maxwell1(node,elem,HB,pde,bdFlag,option); 
    % compute the error
    uI = edgeinterpolate1(pde.exactu,node,edge);
    L2Err(k) = sqrt(abs(real(u-uI)'*eqn.M*real(u-uI)));    
    energyErr(k) = sqrt(abs(real(u-uI)'*eqn.A*real(u-uI)) + L2Err(k)^2);
    L2ErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.M*imag(u-uI)));
    energyErrImag(k) = sqrt(abs(imag(u-uI)'*eqn.A*imag(u-uI)) + L2ErrImag(k)^2);
%     energyErr(k) = getHcurlerror3NE(node,elem,pde.curlu,u);
%     L2Err(k) = getL2error3NE(node,elem,pde.exactu,u);
    N(k) = length(u);
end

%% Plot convergence rates
N = N(1:k); energyErr = energyErr(1:k); L2Err = L2Err(1:k);
energyErrImag = energyErrImag(1:k); L2ErrImag = L2ErrImag(1:k);
figure(1); clf; 
showrate2(N,energyErr,1,'r-+','||Re(u-u_h)||_A',...
          N,L2Err,1,'b-+','||Re(u-u_h)||');
figure(2); clf; 
showrate2(N,energyErrImag,1,'r-+','||Img(u-u_h)||_A',...
          N,L2ErrImag,1,'b-+','||Img(u-u_h)||');