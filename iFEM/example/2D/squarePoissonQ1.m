%% SQUAREPOISSONQ1 Poisson equation in a square domain using bilinear quad element 
%
%   squarePoissonQ1 computes bilinear finite element approximations of the
%   Poisson equation in the unit square on a sequence of quad meshes obtained by
%   uniform refinement. It plots the approximation err vs the number of
%   degree of freedoms.
% 
% See also squarePoisson
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear;

%% Parameters 
maxIt = 4; N = zeros(maxIt,1);
errL2 = zeros(maxIt,1); errH1 = zeros(maxIt,1);

%% Generate an initial mesh 
%node=  [0, 0; 1, 0; 1, 1; 0,1];
node = [0, 0;2 -1; 2.5, 0;2, 1];
elem = [1,2,3,4];
for k = 1:3
    [node,elem] = uniformrefinequad(node,elem);
end


pde = sincosdata;

for k = 1:4
    [uh,Du,eqn,info] = PoissonQ1(node,elem,pde);
    N(k) = size(elem,1);
    
    if N(k) < 2e3 % show mesh and solution for small size
        figure(1);  showsolution(node,elem,uh);    
    end
    errH1(k) = getH1errorQ1(node,elem,pde.Du,uh);  
    errL2(k) = getL2errorQ1(node,elem,pde.exactu,uh);
    
    [node,elem] = uniformrefinequad(node,elem);
end


%% Plot convergence rates
figure(2);
showrate3(N,errH1,1,'-*','||Du-Du_h||',...
          N,errL2,1,'k-+','||u-u_h||');
