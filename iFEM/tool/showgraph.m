function showgraph(node,edge)
%% SHOWGRAPH displays a planar graph
%
%    showgraph(node,edge) displays a planar undirected graph.
%
%   See also showmesh, findedge
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

line([node(edge(:,1),1)'; node(edge(:,2),1)'],...
     [node(edge(:,1),2)'; node(edge(:,2),2)'],...
     'LineWidth',1,'Color',[0.125 0.5 0.125]);
hold on
plot(node(:,1),node(:,2),'k.', 'MarkerSize', 12);
axis equal; axis off