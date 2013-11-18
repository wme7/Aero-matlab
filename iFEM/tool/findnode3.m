function findnode3(node,range,varargin)
%% FINDNODE3 highlights some nodes for 3-D meshes.
%
%    FINDNODE3(node,range) finds all nodes whose indices are in the range
%    array. It plots these nodes using larger black dots and display the indices.
%
%    FINDNODE3(node) finds all nodes.
%   
%    FINDNODE3(node,range,'noindex') skip the display of indices.
%
%    FINDNODE3(node,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the nodes.
%
% Example
%
%     node = [1,1,1; 0,0,0; 1,1,0; 1,0,0];
%     elem = [1,2,3,4];
%     subplot(1,2,1);
%     showmesh3(node,elem);
%     findnode3(node);
%     subplot(1,2,2);
%     showmesh3(node,elem);
%     findnode3(node,1,'noindex','color','r','MarkerSize',24)
%
%   See also findelem3, findnode, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

hold on
if (nargin==1)
    range = (1:size(node,1))'; 
end
if islogical(range) 
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
if (nargin>2) && (strcmp(varargin{1},'coordinate'))
    distance = (node(:,1)-range(1)).^2 + (node(:,2)-range(2)).^2 + ...
               (node(:,3)-range(3)).^2;
    [mind,range] = min(distance); %#ok<ASGLU>
end
h = plot3(node(range,1),node(range,2),node(range,3),'k.','MarkerSize', 20);
if nargin>2
    if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index')
        startidx = 2;
    else
        startidx = 1;
    end
    set(h,varargin{startidx:end});
end
if (nargin<3) || ~(strcmp(varargin{1},'noindex'))
   text(node(range,1)+0.015,node(range,2)+0.015,node(range,3)+0.015,int2str(range), ...
     'FontSize',16,'FontWeight','bold');
end
hold off;