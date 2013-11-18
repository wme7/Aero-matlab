%% Recover function p from P0 to P1
% [rp,nodenew,elemnew] = recoverP02P1(node,elem,p,recoverMethod)
% The recover methods
% 1) 'LS': least square
% 2) 'LA': solve the Laplacian problem, only this method
%           the new node and new elem will be constructed.
%
% Created by Lin Zhong April, 2013. 

%% Node and Elem
clear all; close all;
[node,elem] = squaremesh([0,1,0,1],0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
figure(1)
showmesh(node,elem);
findnode(node);
findelem(node,elem);
display(elem);

NT = size(elem,1);
p = ones(NT,1);
[rp,nodenew,elemnew] = recoverP02P1(node,elem,p,'LA');
figure(2)
showmesh(nodenew,elemnew);
findnode(nodenew);
findelem(nodenew,elemnew);
display(elemnew);

%% LS and LA
% For LS, the value at the nodes are constructed by least square linear
% fitting of the barycenters of the adjacent elemts. 
%
% For LA, we construct new nodes and new elems. The new nodes are the
% vertices of the triangles and the barycenters of the elemenst. The new
% elems are obtained by connecting two barycenters of the interior edges,
% and each end point of that edge. Please check the graph. 
% We get the values of the interior vertices by solving the Laplacian
% problem. To improve the efficiency, we only construct part of the
% stiffness matrix.
% $$ A_{N,N+NT}(rp;p)^t = 0$$
% we solve $$rp = -A_{N,N}\A(N,N+NT)*p$$
% We notice that $A_{N,N}$ is diagonal, so this process should be fast.
% For the boundary nodes, we use least square.
%
% Both methods achieve the same accuracy. 
% Because of the large size for loop, the LS method is slower than LA.
[node,elem] = squaremesh([0,1,0,1],0.05);
NT = size(elem,1);
p = ones(NT,1);
profile on
recoverP02P1(node,elem,p,'LA');
profile viewer


