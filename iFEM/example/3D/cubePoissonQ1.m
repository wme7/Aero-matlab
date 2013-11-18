%% CUBEPOISSONQ1 Poisson equation in a CUBE domain using trilinear cube element 
%
%   cubePoissonQ1 computes trilinear finite element approximations of the
%   Poisson equation in the unit cube on a sequence of hex meshes obtained by
%   uniform refinement. It plots the approximation err vs the number of
%   degree of freedoms.
% 
% See also 
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear all;

%% Parameters 
maxIt = 4; N = zeros(maxIt,1);
errL2 = zeros(maxIt,1); errH1 = zeros(maxIt,1);  erruIuh = zeros(maxIt,1);

h = 1.0/2.0;
cube = [0 1 0 1 0 1];

pde = sincosdata3;
option.solver = 'amg';
% option.isoparametric = true;

for k = 1:maxIt
    [node,elem] = cubehexmesh(cube,h/(2.0^k));
    bdFlag = setboundary3(node,elem,'Dirichlet','abs(z)>eps','Neumann','abs(z)<eps');
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [uh,Du,eqn,info] = Poisson3Q1(node,elem,pde,bdFlag,option);
    [uh,Du,eqn,info] = Poisson3T1(node,elem,pde,bdFlag,option);
    N(k) = size(elem,1);
    uI = pde.exactu(node);
    erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
    errH1(k) = getH1error3Q1(node,elem,pde.Du,uh);  
    errL2(k) = getL2error3Q1(node,elem,pde.exactu,uh);        
end

%% Plot convergence rates
figure;
showrate3(N,errH1,2,'-*','||Du-Du_h||',...
          N,erruIuh,2,'r-+','||Du_I-Du_h||',...
          N,errL2,2,'k-+','||u-u_h||');
      
ts = zeros(k,3); ts = char(ts);
display(' #Dof   ||u-u_h||     ||Du-Du_h||   ||DuI-Du_h||');
display([num2str(N) ts num2str(errL2,'%0.5e') ts num2str(errH1,'%0.5e')...
         ts num2str(erruIuh,'%0.5e')]);