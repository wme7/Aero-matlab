%% 3-D Red Refinement
% 
% It subdivides each tetrahedron in a triangulation into eight
% subtetrahedra of equal volume. 
%
% Reference: J. Bey. Simplicial grid refinement: on Freudenthal's algorithm
% and the optimal number of congruence classes. Numer. Math.. 85(1):1--29,
% 2000. p11 Algorithm: RedRefinement3D.
%  
% S. Zhang. Successive subdivisions of tetrahedra and multigrid methods on
% t etrahedral meshes. Houston J. Math. 21, 541?556, 1995.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Refinement
node = [0,0,0; 1,0,0; 1,1,0; 1,1,1];
elem = [1 2 3 4];
figure(1); subplot(1,2,1);
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.5]);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([34 12]);
findnode3(node);
[node,elem] = uniformrefine3s(node,elem);
figure(1); subplot(1,2,2);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([34 12]);
findnode3(node);
%%
% After cutting the four corner, the remaining octahedron should be divided
% into four tetrahedron by using one of three diagonals. Here follow Bey we
% always use diagonal 6-9. The ordering of sub-tetrahedron is important
% such that recursive application to any initial tetrahedron yields
% elements of at most three congruence classes.

%%
[tempvar,idx] = fixorder3(node,elem);
display(idx);
%%
% Note that the orientation of the 6-th and 8-th children is changed. 

%% Dependence of the Initial Mesh
node = [1,0,0; 1,1,0; 0,0,0; 1,1,1];
elem = [1 2 3 4];
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.5]);
figure(1); subplot(1,2,1);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([34 12]);
findnode3(node);
[node,elem] = uniformrefine3s(node,elem);
figure(1); subplot(1,2,2);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([34 12]);
findnode3(node);
%%
% To have a better mesh quality, one may want to use the diagonal with the
% shortest diagonal (implemented in uniformrefine3l). The current algorithm
% didn't compute the edge length. The mesh quality will depend on the
% ordering of the initial mesh. For example, for the following tet, 6-9 is
% the longest diagonal and the refined mesh is less shape regular although
% still three congruence classes are possible.

%% Uniform mesh of the Cube
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];

figure(1); subplot(1,3,1); 
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.5,0.3]);
showmesh3(node,elem,[],'FaceAlpha',0.25); view([38 10]);
meshquality(node,elem);

[node,elem] = uniformrefine3s(node,elem);
figure(1); subplot(1,3,2);
showmesh3(node,elem,[],'FaceAlpha',0.25); view([38 10]);
meshquality(node,elem);

[node,elem] = uniformrefine3s(node,elem);
figure(1); subplot(1,3,3);
showmesh3(node,elem,[],'FaceAlpha',0.25); view([38 10]);
meshquality(node,elem);
%%
% Starting from a suitable ordering an initial mesh (dividing one cube into
% six tetrahedron), uniformrefine3 (used in cubemesh.m) will produce a
% sequence of uniform mesh of a cube. In the output of mesh quality, only
% one exists which means all tetrahedron are of similar type.

%% Test mesh quality
node = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
elem = [1 2 3 4];
for k = 1:5
    [node,elem] = uniformrefine3s(node,elem);
    meshquality(node,elem);
end
%%
% We test the quality of meshes obtained by uniformrefine3 for a different
% initial mesh. Now the mean of the mesh quality is changing while the
% minimial is bounded below.

%% Test boundary flag
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];
bdFlag = setboundary3(node,elem,'Dirichlet');
for k = 1:3
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    bdFlagnew = setboundary3(node,elem,'Dirichlet');
    display(any(any(bdFlag - bdFlagnew)));
end
%%
% bdFlag obtained by uniformrefine3 is the same as bdFlagnew by finding
% boundary faces of the triangulation.