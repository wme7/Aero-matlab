function [u,A,b,edge] = Poisson3P2(node,elem,bdFace,HB,f,g_D,g_N)
%% POISSON3P2 Poisson equation: P2 quadratic element in 3-D
% u = PoissonP2(node,elem,bdEdge,f,g_D,g_N,d) assembesl the matrix equation
% Au=b for the quadratic finite element discritization of Poisson equation and
% solves it by multigrid solver. The domain is discretized by a
% triangulation represented by |node| and |elem| and the boundary condition
% is given by |bdEdge|. The function handels |f, g_D, g_N| is the data for
% the Poisson equation.
%
% The usage is the same as <a href="matlab:ifemdoc Poisson">Poisson</a>. 
%
% For quadratic elements, middle points of each edge are also degree of
% freedom. See <a href="matlab:ifemdoc P2dof">P2dof</a> for details.
%
% See also Poisson, Poisson3, Poisson3P2 
%
% <a href="matlab:ifemdoc Poisson">ifem Poisson</a>
% <a href="matlab:ifemdoc PoissonP2">ifem PoissonP2</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct Data Structure
[elem2dof,edge] = dof3P2(elem);
N = size(node,1); NT = size(elem,1); Ndof = N + size(edge,1);

%% Compute geometric quantities and gradient of local basis
face = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
v12 = v12(3*NT+1:4*NT,:); 
v13 = v13(3*NT+1:4*NT,:);
v14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = dot(mycross(v12,v13,2),v14,2)/6;
Dlambda1(1:NT,:) = normal(1:NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda2(1:NT,:) = normal(NT+1:2*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda3(1:NT,:) = normal(2*NT+1:3*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda4(1:NT,:) = normal(3*NT+1:4*NT,:)./[6*volume, 6*volume, 6*volume];
% barycentric basis and quadrature points
alpha = 0.58541020; beta = 0.13819660;
lambda1 = [alpha beta beta beta];  % p1 = [a b b b];   w1 = 1/4
lambda2 = [beta alpha beta beta];  % p2 = [b a b b];   w2 = 1/4
lambda3 = [beta beta alpha beta];  % p3 = [b b a b];   w3 = 1/4
lambda4 = [beta beta beta alpha];  % p4 = [b b b a];   w4 = 1/4
weight = 1/4;
% compute element-wise basis
Dphi(:,:,1) = kron(Dlambda1,4*lambda1-1);
Dphi(:,:,2) = kron(Dlambda2,4*lambda2-1);
Dphi(:,:,3) = kron(Dlambda3,4*lambda3-1);
Dphi(:,:,4) = kron(Dlambda4,4*lambda4-1);
Dphi(:,:,5) = 4*(kron(Dlambda1,lambda2) + kron(Dlambda2,lambda1));
Dphi(:,:,6) = 4*(kron(Dlambda1,lambda3) + kron(Dlambda3,lambda1));
Dphi(:,:,7) = 4*(kron(Dlambda1,lambda4) + kron(Dlambda4,lambda1));
Dphi(:,:,8) = 4*(kron(Dlambda2,lambda3) + kron(Dlambda3,lambda2));
Dphi(:,:,9) = 4*(kron(Dlambda2,lambda4) + kron(Dlambda4,lambda2));
Dphi(:,:,10)= 4*(kron(Dlambda3,lambda4) + kron(Dlambda4,lambda3));
clear v12 v13 v14 normal Dlambda1 Dlambda2 Dlambda3 Dlambda4 

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:10
    for j = i:10
        Aij = weight*dot(Dphi(:,:,i),Dphi(:,:,j),2).*volume;
        if (j==i)
            A = A + sparse(double(elem2dof(:,i)),double(elem2dof(:,j)),Aij,Ndof,Ndof);
        else
            A = A + sparse(double([elem2dof(:,i);elem2dof(:,j)]),...
                double([elem2dof(:,j);elem2dof(:,i)]), [Aij; Aij],Ndof,Ndof);                    
        end    
    end
end
clear Aij

%% Assemble right hand side by 5-points quadrature rule
% barycentric coordinate
% p1 = [1/4 1/4 1/4 1/4];   w1 = -4/5
% p2 = [1/2 1/6 1/6 1/6];   w2 = 9/20
% p3 = [1/6 1/2 1/6 1/6];   w3 = 9/20
% p4 = [1/6 1/6 1/2 1/6];   w4 = 9/20
% p5 = [1/6 1/6 1/6 1/2];   w5 = 9/20
lambda1 = [1/4 1/2 1/6 1/6 1/6];
lambda2 = [1/4 1/6 1/2 1/6 1/6];
lambda3 = [1/4 1/6 1/6 1/2 1/6];
lambda4 = [1/4 1/6 1/6 1/6 1/2];
weight = [-4/5 9/20 9/20 9/20 9/20];
% x-y-z coordinate
nQuad = length(lambda1);
pxy = zeros(NT,3,nQuad);
for p = 1:nQuad
    pxy(:,:,p) = lambda1(p)*node(elem(:,1),:) ...
               + lambda2(p)*node(elem(:,2),:) ...
               + lambda3(p)*node(elem(:,3),:) ...
               + lambda4(p)*node(elem(:,4),:);
end
%----------- Assembing right hand side by 5-point rule --------------------
phi(1,:) = lambda1.*(2*lambda1-1);
phi(2,:) = lambda2.*(2*lambda2-1);
phi(3,:) = lambda3.*(2*lambda3-1);
phi(4,:) = lambda4.*(2*lambda4-1);
phi(5,:) = 4*lambda1.*lambda2;
phi(6,:) = 4*lambda1.*lambda3;
phi(7,:) = 4*lambda1.*lambda4;
phi(8,:) = 4*lambda2.*lambda3;
phi(9,:) = 4*lambda2.*lambda4;
phi(10,:)= 4*lambda3.*lambda4;
bt = zeros(NT,10);
for j = 1:10
    for p = 1:nQuad
        bt(:,j) = bt(:,j) + f(pxy(:,:,p)).*phi(j,p)*weight(p);
    end
    bt(:,j) = bt(:,j).*volume;
end
b = accumarray(elem2dof(:),bt(:),[Ndof 1]);

%% Boundary Conditions
% Find boundary edges and nodes
if ( isempty(bdFace) && ~isempty(g_D))
    bdFace = setboundary3(elem,1);
end
Dirichlet = [elem2dof((bdFace(:,1)==1),[2 4 3 10 8 9]); ...
             elem2dof((bdFace(:,2)==1),[1 3 4 10 7 6]); ...
             elem2dof((bdFace(:,3)==1),[1 4 2 9 5 7]);...
             elem2dof((bdFace(:,4)==1),[1 2 3 8 6 5])];
isBdDof = false(Ndof,1); 
isBdDof(Dirichlet) = true;
bdDof = find(isBdDof);
% freeDof = find(~isBdDof);

if ( isempty(bdFace) && ~isempty(g_N))
    bdFace = setboundary3(elem,2);
end
isBdface = reshape(bdFace,4*NT,1);
Neumann = face((isBdface == 2),:); 
%-------------------- Dirichlet boundary conditions------------------------
u = zeros(Ndof,1);
if ~isempty(g_D)
    idx = (bdDof > N);
    u(bdDof(~idx)) = g_D(node(bdDof(~idx),:));
    bdEdgeIx = bdDof(idx)-N;
    xyzbdEdgeDof = (node(edge(bdEdgeIx,1),:) + node(edge(bdEdgeIx,2),:))/2;
    u(bdDof(idx)) = g_D(xyzbdEdgeDof);    
    b = b - A*u;
    b(bdDof) = u(bdDof);
end
% % %-------------------- Neumann boundary conditions -----------------------
if (~isempty(Neumann) && ~isempty(g_N))
    v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
    v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
    area = 0.5*sqrt(sum(mycross(v12,v13,2).^2,2));
    centroid = (node(Neumann(:,1),:)+node(Neumann(:,2),:) ...
               +node(Neumann(:,3),:)+node(Neumann(:,4),:))/4;
    gt = area.*g_N(centroid)/4;
    b = b + accumarray([Neumann(:),ones(3*size(Neumann,1),1)], ... 
                       repmat(gt,3,1),[N,1]);
    if  isempty(g_D) % compatible condition
        b = b - mean(b);
    end
end
if  (isempty(g_D) && isempty(g_N)) % pure Neumann boundary condition
        b = b - mean(b);
%         freeDof = 2:Ndof; % pure Neumann boundary condition has a kernel
end

%% Solve the linear system of algebraic equations
bdidx = zeros(Ndof,1); 
bdidx(bdDof) = 1;
Tbd = sparse(1:Ndof,1:Ndof,bdidx,Ndof,Ndof);
T = sparse(1:Ndof,1:Ndof,1-bdidx,Ndof,Ndof);
AD = T*A*T + Tbd;
%     tic; u = AD\b; display('direct solver'); toc
u = mg(AD,b,elem,[],HB,edge);
