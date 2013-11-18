function examplePoisson3CR

clc; clear all; close all;
maxIt = 4; N = zeros(maxIt,1);
% erruIuh = zeros(maxIt,1);
err = zeros(maxIt,1);

node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];  % nodes
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7]; % elements

pde = sincosdata3;
for k = 1:maxIt
    [u,Duh,eqn] = Poisson3CRrevise(node,elem,pde);
%     A = eqn.A;
%     face = eqn.face;
    N(k) = size(elem,1);
%     uI = pde.exactu((node(face(:,1),:)+node(face(:,2),:)+ ...
%                      node(face(:,3),:))/3);
%     erruIuh(k) = sqrt((u-uI)'*A*(u-uI));
    err(k) = getH1error3(node,elem,pde.Du,Duh);
    [node,elem] = uniformrefine3(node,elem);
end

r = showrate(N,err,2,'-*');
legend('||DuI-Duh||',['N^{' num2str(r) '}'], 'LOCATION','Best');