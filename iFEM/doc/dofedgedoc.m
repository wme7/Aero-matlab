%% Dofedge in Two Dimensions
% 
% We describe general idea of the data structures generated in subroutine 
% dofedge for two dimensional triangular grids, suppose $N,NE,NT$ represents
% the number of vertices, edges and element respectively. 

%% Input
% |elem| is $NT\times 3$ matrix with $NT$ be the number of
% the triangles.  |elem(i,:)| contain the global index of the vertex of the
% |i-th| triangle.

%% Output
% * |elem2dof|, a $NT\times 3$ matrix, is the elementwise pointer from elem 
% to dof indices. |elem2dof(i,:)| stores the global indices of edges for triangle |i|.
% * |dofsign|, $NT\times 3$ matrix with element 1 or -1, records the consistency of the local and global
% edge orientation. value 1 means the global and local direction is the
% same while -1 means the oposite.
% * |edge| is the edge matrix and edge(:,1)<edge(:,2). The orientation of
% edge is from the node with smaller global index to bigger one. In each
% triangle, the local edge orientated as [2 3; 3 1; 1 2].
% * |edgeSign|, $NE\times 1$, indicates the consistence between the global orientation of an
% edge and the local orientation of the first element containing this
% edge.It mainly used to deal with the boudary condition. When the global
% direction of the boudary edge is the same with the domain direction,we set the value
% 1, otherwise, -1.

%% Example
%Generate edge and elem2dof
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 5,6,8; 7,8,6];
node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1;0,-1;0,0];
T = auxstructure(elem);
elem2dof = T.elem2edge;
edge = T.edge;
edge2elem = T.edge2elem;

showmesh(node,elem);
findnode(node,[],'index');
findelem(node,elem,[],'index');
findedge(node,edge);

display(elem);
display(elem2dof);
display(edge);


% Consistency of oritentation of edges
NT = size(elem,1); NE = size(edge,1);
edgeSign = -ones(NE,1); 
isConsistentEdge = (elem(edge2elem(:,1)+NT*mod(edge2elem(:,3),3))== edge(:,1)); 
edgeSign(isConsistentEdge) = 1;

% The sign of dof
isbdEdge = (edge2elem(:,1) == edge2elem(:,2));
dofSign = zeros(NT,3,'int8');
dofSign(edge2elem(:,1) + NT*(edge2elem(:,3)-1)) = edgeSign;
dofSign(edge2elem(~isbdEdge,2)+NT*(edge2elem(~isbdEdge,4)-1)) = -edgeSign(~isbdEdge);

display(edgeSign);
display(dofSign);

%% Detail explaination of the code
% Matrix |isConsistentEdge| indicates the consistent direction of the edges. We
% treat structure |elem| here as the $3NT\times 1$ column vector.
% |t = edge2elem(e,1)| means that |t| is the first triangle containing |e|.
%
% |k = edge2elem(e,3)| means the |k-th| edge of |t| is |e|.
%
% |elem(edge2elem(:,1) + NT*mod(edge2elem(:,3),3))| is the first node of
% the |k|-th edge of |t|. 
%
% |elem(edge2elem(:,1)+NT*mod(edge2elem(:,3),3))==edge(:,1)| means when the 
% first node of the |k-th| edge of |t| is the same as
% that in the global index, then, the local and global directions of the edge
% is consistent. |edgeSign| equals 1 if the global direction of edge is 
% consistent with the local direction of edge in its first element, equals -1 otherwise.

