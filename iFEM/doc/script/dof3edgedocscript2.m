node = [1,0,0; 0,1,0; 0,0,0; 0,0,1];
elem = [2 4 3 1];
localEdge = [elem(:,1) elem(:,2); elem(:,1) elem(:,3); elem(:,1) elem(:,4); ...
             elem(:,2) elem(:,3); elem(:,2) elem(:,4); elem(:,3) elem(:,4)];
edge = sort(localEdge,2);
elem2dof = [1 2 3 4 5 6];
dofSign = [1 1 -1 -1 -1 -1];     
figure(1); clf;
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.6,0.4]);
subplot(1,2,1)
showmesh3(node,elem);
view(-14,12);
tempnode = node(elem(:),:);
findnode3(tempnode);
findedge3(node,localEdge);
midEdge = (node(localEdge(:,2),:) + node(localEdge(:,1),:))/2;
edgeVec = node(localEdge(:,2),:) - node(localEdge(:,1),:);
h = quiver3(midEdge(:,1),midEdge(:,2),midEdge(:,3), ...
            edgeVec(:,1),edgeVec(:,2),edgeVec(:,3));
set(h,'Linewidth',2)
subplot(1,2,2)
showmesh3(node,elem);
view(-14,12);
findnode3(node);
findedge3(node,edge);
midEdge = (node(edge(:,2),:) + node(edge(:,1),:))/2;
edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
h = quiver3(midEdge(:,1),midEdge(:,2),midEdge(:,3), ...
            edgeVec(:,1),edgeVec(:,2),edgeVec(:,3));
set(h,'Linewidth',2)