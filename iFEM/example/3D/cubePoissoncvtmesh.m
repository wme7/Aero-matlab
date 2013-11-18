function [N,errL2,errH1,erruIuh] = cubePoissoncvtmeshnew
% CUBEPOISSON solves Poisson equation in a cube.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all; clear all;
%% Parameters
maxIt = 4; N = zeros(maxIt,1); 
errL2 = zeros(maxIt,1); errH1 = zeros(maxIt,1); erruIuh = zeros(maxIt,1);

%% Get the data of the pde
pde = sincosdata3;
% pde = polydata3;

%% Finite Element Method        
for k = 1:maxIt
    [node,elem] = cvtuniformmesh(2^k);
    % solve the equation
    option.solver = 'direct';
    [u,Du,A,eqn] = Poisson3(node,elem,pde,[],[],option); 
    N(k) = size(node,1);
    % compute error
    errL2(k) = getL2error3(node,elem,pde.exactu,u);
    errH1(k) = getH1error3(node,elem,pde.Du,Du);
    uI = pde.exactu(node);  % nodal interpolation
    erruIuh(k) = sqrt((u-uI)'*eqn.AD*(u-uI));
end

%% Plot convergence rates
figure(2);
r1 = showrate(N,errH1,2,'-*');
hold on;
r2 = showrate(N,errL2,2,'k-+');
r3 = showrate(N,erruIuh,2,'m-+');
legend('||Du-Du_h||',['N^{' num2str(r1) '}'], ...
       '||u-u_h||',['N^{' num2str(r2) '}'], ...
       '||DuI-Du_h||',['N^{' num2str(r3) '}'], 'LOCATION','Best');
%%
% The error in H1 and L2 norm converges at optimal rate. But no
% superconvergence on criss-cross grids.
end