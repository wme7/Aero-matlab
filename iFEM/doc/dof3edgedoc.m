%% Data Structure: Lowest Order Edge Element
% We explain degree of freedoms for the lowest order edge element on
% tetrahedrons. The required data structure can be constructured by
%
%    [elem2dof,dofSign,edge] = dof3edge(elem)
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Local Labeling of DOFs
dof3edgedocscript1
%%
% The six dofs associated to edges in a tetrahedron is sorted in the
% ordering [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]. Here [1 2 3 4] are local indices
% of vertices.
%
% For a triangulation consisting of several tetrahedrons, all local edges
% are collected and the duplication is eliminated to form edges of the
% triangulation. Therefore the global indices of edges is different with
% the local one. elem2dof(:,1:6) records the mapping from local to global
% indices.
display(elem2dof)

%% Orientation of Edges
dof3edgedocscript2
%%
% The edge [i,j] is orientated in the direction such that i<j, where i,j
% are *global* indices. The edges formed by local indices may not be
% consistent with this orientation. Therefore dofSign is used to indicate
% the consistency.
%
% The nodal indices in the left figure is local while that in the right is
% the global one. The local direction and global direction of edges is also
% indicated by different edge vectors.
display(elem2dof)
display(dofSign)

%% An Example
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
[elem2dof,dofSign,edge] = dof3edge(elem);
figure(1); clf;
showmesh3(node,elem,[],'FaceAlpha',0.25);
view(-30,10);
findedge3(node,edge);
findelem3(node,elem,[1 3]');
display(elem2dof);
display(dofSign);
