%% COARSEN in THREE DIMENSIONS


%% Output
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
bdFlag = setboundary3(node,elem,'Dirichlet');
[node,elem,bdFlag,HB] = bisect3(node,elem,1,bdFlag,HB);
display('Original mesh');
display(elem);
display(bdFlag);
display(HB);
%%
% The node 9 is the new node added by the bisection and is the last index
% of each element by the bisect rule. HB(9,2:3) are two parent nodes of 9
% and HB(9,4) is the generation. All boundary faces are opposite to node 9
% in each element.

figure(2); subplot(1,2,1);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([210 8]);
findnode3(node);
[node,elem,HB,bdFlag,tree,indexMap] = coarsen3(node,elem,'all',HB,bdFlag);
figure(2); subplot(1,2,2);
showmesh3(node,elem,[],'FaceAlpha',0.15); view([210 8]);
findnode3(node);
display(elem);
display(bdFlag);
display(HB);
display(tree);
display(indexMap);

%% Uniform Coarsen
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],1);
[node,elem,~,HB] = uniformbisect3(node,elem,[],HB);
figure(1); subplot(1,3,1); 
set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.65,0.4]);
showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
[node,elem,HB] = coarsen3(node,elem,'all',HB);
figure(1); subplot(1,3,2);
showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
[node,elem,HB] = coarsen3(node,elem,'all',HB);
figure(1); subplot(1,3,3);
showmesh3(node,elem,[],'FaceAlpha',0.35); view([210 8]);
%%
% coarsen3 with all elements marked will remove only half of the elements
% while uniformbisect3 will bisect one element three times.