%% Local Labeling of DOFs
node = [1,0,0; 0,1,0; 0,0,0; 0,0,1];
elem = [1 2 3 4];
localEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
edge = zeros(20,2);
edge([1 12 5 20 11 4],:) = localEdge;
elem2dof = [1 12 5 20 11 4];
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.6,0.4]);
subplot(1,2,1)
showmesh3(node,elem);
view(-14,12);
findnode3(node);
findedge3(node,localEdge);

%% Face and face to local DOF mapping
face1 = [2 3 4];
f1dof = [6 5 4];
face2 = [1 4 3];
f2dof = [6 2 3];
face3 = [1 2 4];
f3dof = [5 3 1];
face4 = [1 3 2];
f4dof = [4 1 2];

%% 