
%% Problem set up
% Lshape domain
[node,elem] = squaremesh([-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
pde = Lshapedata; 
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 2;

%% Finite Element Method
femPoisson(node,elem,pde,bdFlag,option);

%% Adaptive Finite Element Method
afemPoisson(node,elem,pde,bdFlag,option);