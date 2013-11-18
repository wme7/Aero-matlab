%% Coloring Vertices of a Graph
%
% nodeSet = coloring(G,nc) color the vertices of the graph G by nc colors.
% Vertices in the same color are disjointed except the last one. The input
% graph G is represented by the sparse pattern of a symmetric sparse
% matrix. The values of G is not used. The output nodeSet is a cell
% array and nodeSet{k} is the indices of vertices in the k-th color. By
% default, nc = 6 if not specified by the input.
%
% nodeSet = coloring(G,nc,qp) chose the vertices with a larger value qp
% first. For example, qp can be the degree of vertices. Then nodeSet{1}
% will pick up nodes with a local max degree. If no qp, random values are
% used.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear all; close all
%% Generate a test graph
load lakemesh
N = size(node,1); 
NT = size(elem,1);
t2p = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
G = t2p'*t2p;  % connectness of nodes 

%% Set up
N = size(G,1);
N0 = min(floor(N/100),20);   % number of the coarest nodes
if ~exist('nc','var'), nc = 6; end
nodeSet = cell(nc,1);
if ~exist('qp','var')
    qp = rand(N,1); 
else
    qp = (1 + 0.1*rand(N,1))*qp;    
end

%% Coloring
colorstring = ['r' 'y' 'b' 'm' 'k' 'g'];
for k = 1:nc-1
    isC = false(N,1);       % C: coarse nodes
    isF = false(N,1);       % F: fine nodes
    isU = (qp>0);           % U: undecided nodes
    if ~any(isU)            % no more vertices left
        break
    end
    for m = 1:round(log10(N))
        % Mark all available nodes
        isS = isU;          % selected all undecided vertices
        S = find(isS); 

        % Find marked nodes with local maximum quality
        [i,j] = find(triu(G(S,S),1));   % i,j and i<j: edges of subgraph S
        idx = (qp(S(i)) >= qp(S(j)));   % compare quality of connected vertices
        isS(S(j(idx))) = false;         % remove vertices with smaller qp
        isS(S(i(~idx))) = false;        % if qualities are equal, keep the nodes with smaller index
        isC(isS) = true;                % set selected nodes as coarse nodes

        % Remove coarse nodes and neighboring nodes from undecided set
        [i,j] = find(G(:,isC)); %#ok<*NASGU>
        isF(i) = true;                  % neighbor of C nodes are F nodes
        isU(isF | isC) = false;         % remove current C and F from U
        if sum(isU) <= N0               
            break;    % only break the inner loop
        end
    end
    nodeSet{k} = find(isC);         
    qp(nodeSet{k}) = 0;    
end
nodeSet{k+1} = find(qp>0);  % all left vertices are added to the last color

%% Show coloring of nodes
set(gcf,'Units','normal'); 
set(gcf,'Position',[0,0,0.6,0.6]);
showcoloring(G,node,nodeSet);

%%
% Since all left vertices are added to the last color, vertices in the last
% one (in green) are not disjoint. But usually it is of small size and cost but robust
% algorithms can be further applied to coloring it. In the application of
% iterative methods, subproblems restricted to the last color can be solved
% exactly since the dimension is small.