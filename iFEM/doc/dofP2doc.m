%% Data Structure: P2 Quadratic Element
%
% We explain degree of freedoms for quadratic element on triangles. There
% are two types of dofs: vertex type and edge type. Given a mesh, the
% required data structure can be constructured by
%
%   [elem2dof,edge,bdDof] = dofP2(elem)

help dofP2

%% Local indexing of DOFs
node = [1,0; 1,1; 0,0];
elem = [1 2 3];
edge = [2 3; 1 3; 1 2];
% elem2dof = 1:6;
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.5,0.3]);
subplot(1,2,1)
showmesh(node,elem);
findnode(node);
findedgedof(node,edge);
subplot(1,2,2)
showmesh(node,elem);
findnode(node);
findedge(node,edge);
%%
% The six dofs in a triangle is displayed in the left. The first three are
% associated to the vertices of the triangle and the second three to the
% middle points of three edges. The edges are indexed such that the
% i-th edge is opposite to the i-th vertex for i=1,2,3.

%% Global indexing of DOFs
node = [0,0; 1,0; 1,1; 0,1];
elem = [2,3,1; 4,1,3];      
[node,elem] = uniformbisect(node,elem);
figure(2); clf;
showmesh(node,elem);
findnode(node);
findelem(node,elem);
[elem2dof,bdDof,edge] = dofP2(elem);
findedgedof(node,edge);
%%
% The matrix elem2dof is the local to global index mapping of dofs. The
% first 3 columns, elem2dof(:,1:3), is the global indices of dofs
% associated to vertexes and thus elem2dof(:,1:3) = elem. The columns
% elem2dof(:,4:6) point to indices of dofs on edges. The matrix bdDof
% contains all dofs on the boundary of the mesh.
display(elem2dof);
display(bdDof);

%% Local bases of P2
node = [0,0; 1,0; 0,1];
elem = [1 2 3];
edge = [2 3; 1 3; 1 2];
% elem2dof = 1:6;
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.3,0.3]);
showmesh(node,elem); axis on
findnode(node);
findedgedof(node,edge);
%%
% The six Lagrange-type bases functions are denoted by
% $\phi_i, i=1:6$, i.e. $\phi_i(x_j)=\delta _{ij},i,j=1:6$. In
% barycentric coordinates, they are 
%
% $$ \phi_1 = \lambda_1(2\lambda_1 -1),\quad \nabla \phi_1 = \nabla \lambda_1(4\lambda_1-1),$$
%
% $$ \phi_2 = \lambda_2(2\lambda_2 -1),\quad  \nabla \phi_2 = \nabla \lambda_2(4\lambda_2-1),$$ 
%
% $$ \phi_3 = \lambda_3(2\lambda_3 -1),\quad  \nabla \phi_3 = \nabla \lambda_3(4\lambda_3-1),$$ 
%
% $$ \phi_4 = 4\lambda_2\lambda_3,\quad  \nabla\phi_4 = 4\left (\lambda_2\nabla \lambda_3 + \lambda_3\nabla \lambda_2\right )\; ,$$ 
%
% $$ \phi_5 = 4\lambda _3\lambda_1,\quad  \nabla\phi_5= 4\left (\lambda_3\nabla \lambda_1 + \lambda_1\nabla \lambda_3\right )\; ,$$ 
%
% $$ \phi_6 = 4\lambda _1\lambda_2,\quad  \nabla\phi_6=4\left (\lambda_1\nabla
% \lambda_2 + \lambda_2\nabla\lambda_1\right )\; .$$
% 
% When transfer to the reference triangle formed by $(0,0),(1,0),(0,1)$,
% the local bases in x-y coordinate can be obtained by substituting 
% 
% $$\lambda _1 = x, \quad \lambda _2 = y, \quad \lambda _3 = 1-x-y.$$ 
