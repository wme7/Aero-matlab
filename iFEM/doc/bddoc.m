%% Data Structure: Boundary Conditions
%
% We use |bdFlag(1:NT,1:3)| to record the type of three edges of each
% triangle. Similarly in 3-D, we use |bdFlag(1:NT,1:4)| to record the type
% of four faces of each tetrahedron. The value is the type of boundary
% condition.
%
% * 0: non-boundary, i.e., an interior edge or face.
% * 1: first type, i.e., a Dirichlet boundary edge or face. 
% * 2: second type, i.e., a Neumann boundary edge or face. 
% * 3: third type, i.e., a Robin boundary edge or face.

%% Local labeling of edges and faces
% We label three edges of a triangle such that |bdFlag(t,i)| is the edge
% opposite to the i-th vertex. Similarly |bdFlag(t,i)| is the face opposite
% to the i-th vertex.

node = [1,0; 1,1; 0,0];
elem = [1 2 3];
edge = [2 3; 1 3; 1 2];
showmesh(node,elem);
findnode(node);
findedge(node,edge);

%% Set up boundary conditions
%
% The function |setboundary| is to set up the bdFlag matrix for a 2-D
% triangulation and |setboundary3| for a 3-D triangulation. 
%
help setboundary

%% 
% Note that if the i-th edge of t is on the boundary but |bdFlag(t,i)=0|,
% it is equivalent to use homogenous Neumann boundary condition (zero
% flux).


%% Example: Crack Domain
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
elem = label(node,elem);                        % label the mesh
figure;
showmesh(node,elem);                            % plot mesh
findelem(node,elem);                            % plot element indices
findnode(node,2:6);                             % plot node indices
text(node(6,1),node(6,2)+0.075,int2str(1),'FontSize',16,'FontWeight','bold');
hold on;
plot([node(1,1),node(5,1)], [node(1,2),node(5,2)],'r-', 'LineWidth',3);
bdFlag = setboundary(node,elem,'Dirichlet');                % Dirichlet boundary condition
display(elem)
display(bdFlag)
%% 
% node 1 and node 6 are the same point (1,0)

%%
bdFlag = setboundary(node,elem,'Dirichlet','abs(x) + abs(y) == 1','Neumann','y==0');
display(bdFlag)

%% Example: Prism Domain
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
elem = [1,2,3,7; 1,6,2,7; 1,5,6,7];
elem = label3(node,elem);
figure;
showmesh3(node,elem);
view([-53,8]);
findnode3(node,[1 2 3 5 6 7]);
findelem3(node,elem);
bdFlag = setboundary3(node,elem,'Dirichlet','(z==1) | (z==-1)');
display(elem)
display(bdFlag)

%%
% The top and bottom of the prism is set as Dirichlet boundary condition
% and other faces are zero flux boundary condition.

%% Remark
% It would save storage if we record boundary edges or faces only. The
% current data structure is convenient for the local refinement and
% coarsening since the boundary can be easily update along with the change
% of elements. The matrix |bdFlag| is sparse but we use a dense matrix to
% store it. We do not save |bdFlag| as a sparse matrix since updating
% sparse matrix is time consuming. We set up the type of |bdFlag| or
% |bdFlag| to |uint8| to minimize the waste of spaces.

