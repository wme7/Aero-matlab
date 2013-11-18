function [node,elem,u] = ellipsoidPoisson

maxIt = 5;
err = zeros(maxIt,1);
% ----------------------- generate initial mesh ---------------
node = [1,0,0; 0,1,0; -1,0,0; 0,-1,0; 0,0,1; 0,0,-1];
elem = [6,1,2; 6,2,3; 6,3,4; 6,4,1; 5,1,4; 5,3,4; 5,3,2; 5,2,1];
% node = [1,0,0; 0,1,0; -1,0,0; 0,-1,0; 0,0,1];
% elem = [5,1,4; 5,3,4; 5,3,2; 5,2,1];
showmesh(node,elem)
for i=1:3
    [node,elem] = uniformrefine(node,elem);
    %[node,elem] = uniformbisect(node,elem);
end
% project the mesh to the ellipsoid
r = sqrt(node(:,1).^2+node(:,2).^2+node(:,3).^2);
node = node./[r r r];
node(:,1) = 3*node(:,1); node(:,2) = 2*node(:,2);
showmesh(node,elem);

for k=1:maxIt
%    [u,A] = surfacePoisson(node,elem,[],@f,@g_D,[]);
    [u,A] = surfacePoisson(node,elem,[],@f,[],[]);
    uI = exactu(node);
    err(k) = sqrt(u-uI)'*A*(u-uI);
    [node,elem] = uniformbisect(node,elem);
    r = sqrt(node(:,1).^2+node(:,2).^2+node(:,3).^2);
    node = node./[r r r];
    node(:,1) = 3*node(:,1); node(:,2) = 2*node(:,2);
end
plot(err)
end

%---------------- Data of PDE-----------------------------
function z = f(p) % load data (right hand side function)
z = 2*p(:,1);
end
%--------------------------------------------------------------------------
function z = g_D(p) % Dirichlet boundary condition
z = exactu(p);
end
%--------------------------------------------------------------------------
function z = exactu(p)   % exact solution
z = p(:,1);
end
%--------------------------------------------------------------------------