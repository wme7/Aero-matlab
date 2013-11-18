function findedge(node,edge,range,varargin)
%% FINDEDGE highlights edges
%
%    FINDEDGE(node,edge,range) finds all elements whose indices are in the
%    range array by displaying these elements in yellow.
%
%    FINDEDGE(node,edge) finds all elements.
%   
%    FINDEDGE(node,edge,'noindex') skip the display of indices.
%
%    FINDEDGE(node,edge,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the elements.
%    
% Example:
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];
%     subplot(1,2,1);
%     showmesh(node,elem);
%     T = auxstructure(elem);
%     findedge(node,T.edge,5,'index','color','r','MarkerSize',24);
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findedge(node,T.edge);
%
%   See also findelem3, findnode3, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% set up range
if (nargin<=2) || (isempty(range))
    range = (1:size(edge,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
% draw edges in red
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
hold on
if length(range) < size(edge,1)
    h = line([node(edge(range,1),1)'; node(edge(range,2),1)'],...
             [node(edge(range,1),2)'; node(edge(range,2),2)'],...
             'LineWidth',2,'Color','r');
 else
    h = plot(midEdge(:,1),midEdge(:,2),'r.','MarkerSize', 18);
end
% else
%     h = plot(midEdge(:,1),midEdge(:,2),'r.','MarkerSize', 18);
% end
if nargin > 3
    if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index') && length(varargin)>1
        startidx = 2;
    else
        startidx = 1;
    end
    set(h,varargin{startidx:end});
end
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    text(midEdge(:,1)+0.025,midEdge(:,2)-0.025,int2str(range), ...
         'FontSize',12,'FontWeight','bold','Color','k');
end