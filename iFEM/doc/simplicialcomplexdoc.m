%% Simplicial Complex in 3D
%
% We describe general idea of the data structure generated for three
% dimensional triangular grids, suppose $N,NE,NF,NT$ represents the number
% of vertices, edges, faces and elements resprectively. 

%% Construct Data Structure
% elem is $NT\times 4$ matrix with $NT$ be the number of the tetrahedrons.
% elem(i,:) contain the global index of the vertex of the i-th tetrahedron.
% we sort the vertices of elem such
% that:elem(i,1)<elem(i,2)<elem(i,3)<elem(i,4). In this way, the volume is
% not always positive, because we sacrifice the righthand orientation of tet.
% We also sort vertices of face such that:face(j,1)<face(j,2)<face(j,3),and
% the vertices of edge such that: edge(k,1)<edge(k,2).
%
% The local ordering of face in a tetrahedron is locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3].
% The local ordering of edge in a face is locEdge = [2 3;1 3;1 2].
% Now the ordering of local indices is always consistent with that of global
% indices.
% 
% In the following example, elem2face, a $NT\times 4$ matrix, is the
% elementwise pointer from elem to face indicies. elem(i,:) stores the
% global indicies of faces for tetrahedron i. elem2edge, a $NT\times 6$ matrix, is the
% elementwise pointer from elem to edge indicies. elem(j,:) stores the
% global indicies of edges for tetrahedron j. 

%% Example 
% Generate 

elem = [1 4 5 8;1 4 5 7];
node = [1,0,0;1,1,1;1,-1,-1; 0,1,0; 0,0,0;1,1,-1; 0,0,1;-1,-1,-1];
[elem2edge,~,edge] = dof3edge(elem);
[elem2face,~,face] = dof3RT0(elem);
figure(1);clf;
showmesh3(node,elem);
view(-20,20);
findnode3(node,[1,4,5,7,8]);
findelem3(node,elem);
display(edge); display(elem2edge);
display(face); display(elem2face);

%% The directions. 
%% Direction of face.
% While we sort the elem, face and edge in the increasing order of
% vertices. We determine the direction of a face [i j k] by the righthand
% rule. The direction of face is not always outword normal direction of  
% a certain tetrahedron. We denote the direction as +1 if the direction
% of a face is the same with the normal direction in a certain elem, and -1
% otherwise. Then for locface of a tetrahedron, the directions should be
% sign(volume)*[+1 -1 +1 -1].
%% Direction of edge.
% Similarly, we determine the direction of a edge [i j] is from the node with 
% smaller global index to bigger one(that is from i to j). The direction of
% edge is not always consistent with the counter clockwise direction in a
% face. When the edge direction is the same with the local counter
% clockwise direction, we set the direction as +1, otherwise -1. Then for
% locedge of a face, the directions should be [+1 -1 +1].
%% Direction of the basis.
% According to our sorted elem and face, the RT0 basis on face [i j k].
% Defined as $\phi_{i,j,k} = 2(\lambda_i \bigtriangledown \lambda_j \times
% \bigtriangledown \lambda_k+
% \lambda_j \bigtriangledown \lambda_k \times \bigtriangledown
% \lambda_i+\lambda_k \bigtriangledown \lambda_i \times \bigtriangledown
% \lambda_j)$ The direction is always consistent with the face direction.
%
%% Inward normal direction.
% Dlambda(t,:,1) is the gradient of lambda of the 1st index of the
% t-th tet. The direction of Dlambda is always inward normal direction. 
[Dlambda, ~] = gradbasis3(node,elem(2,:)); 
display(Dlambda);



