% function [N, erruIuh, r1] = squareelasticity
%  
close all;clc; clear;
%% Parameters 
maxIt = 4; N = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);


%% Generate an initial mesh 
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
para.mu = 1;
para.lambda = 1;
pde = elasticitydata(para);
for i = 1:6
   [node,elem] = uniformrefine(node,elem); 
end

bdEdge = setboundary(node,elem,'Dirichlet');
for i = 1:maxIt

    [u,A,Mf] = elasticity(node,elem,pde,bdEdge);
    uI = pde.exactu(node);
    uI = uI(:);
    erruIuh(i) = sqrt((u - uI)'*A*(u - uI));
    N(i) = size(node,1);
    [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
end

figure(2);
r1 = showrate(N,erruIuh,2,'k-*');