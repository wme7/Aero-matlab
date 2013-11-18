%% Data Structure: CR Nonconforming P1 Element
%
% We explain degree of freedoms for CR nonconforming P1 element on
% triangles and tetrahedrons. The dofs are associated to edges (2-D) and
% faces (3-D). Given a mesh, the required data structure can be
% constructured by
%
%   [elem2edge,edge] = dofedge(elem);
%   [elem2face,face] = dof3face(elem);

%% Local indexing of DOFs
node = [1,0; 1,1; 0,0];
elem = [1 2 3];
edge = [2 3; 1 3; 1 2];
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.5,0.3]);
subplot(1,2,1)
showmesh(node,elem);
findnode(node);
findedge(node,edge);
node = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
elem = [1 2 3 4];
face = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
subplot(1,2,2)
showmesh3(node,elem); 
view([-26 10]);
findnode3(node);
findelem(node,face);

%%
% The three dofs associated to edges in a triangle is displayed in the left
% and the four dofs associated to faces of a tetrahedron is in the right.
% The dofs are indexed such that the i-th dof is opposite to the i-th
% vertex.

%% Local bases of CR element
%
% The d+1 Lagrange-type bases functions are denoted by
% $\phi_i, i=1:d+1$, i.e. $\phi_i(m_j)=\delta _{ij},i,j=1:d+1$. In
% barycentric coordinates, they are:
%
% * 2-D:  $\phi_i = 1- 2\lambda_i,\quad \nabla \phi_i = -2\nabla \lambda_i,$
% * 3-D:  $\phi_i = 1- 3\lambda_i,\quad \nabla \phi_i = -3\nabla \lambda_i.$

%% Global indexing of DOFs
node = [0,0; 1,0; 1,1; 0,1];
elem = [2,3,1; 4,1,3];      
[node,elem] = uniformbisect(node,elem);
figure(2); clf;
showmesh(node,elem);
findnode(node);
findelem(node,elem);
[elem2edge,edge] = dofedge(elem);
findedge(node,edge);
%%
% The matrix elem2edge is the local to global index mapping of edges.
display(elem2edge);