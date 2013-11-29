% (1) compute connectivity arrays for triangles
function [ElmtToElmt, ElmtToEdge] = umCONNECT(ElmtToNode)
Nedges = 3;

% (2) find nber of elements, maximum node number, total local edges
K = size(ElmtToNode,1);
Nnodes = max(max(ElmtToNode));
TotalEdges = Nedges*K;

% (3) list of local Edge to local Node connections
vnum = [[1,2];[2,3];[1,3]];

% (4) build global Edge to node  array
EdgeToNode = spalloc(TotalEdges, Nnodes, 2*TotalEdges);
sk = 1;
for k=1:K
    for Edge=1:Nedges
        EdgeToNode( sk, ElmtToNode(k, vnum(Edge,:))) = 1;
        sk = sk+1;
    end
end

% (5) build global Edge to global Edge  array
EdgeToEdge = EdgeToNode*transpose(EdgeToNode);

% (6) remove self connections
EdgeToEdge = EdgeToEdge - 2*speye(TotalEdges);

% (7) find complete Edge to Edge connections
[Edges1, Edges2] = find(EdgeToEdge==2);

% (8) convert Edge global number to element and local Edge numbers
element1 = floor( (Edges1-1)/Nedges )  + 1;
Edge1    =   mod( (Edges1-1), Nedges ) + 1;
element2 = floor( (Edges2-1)/Nedges )  + 1;
Edge2    =   mod( (Edges2-1), Nedges ) + 1;
% (9) rearrange into K x Nedges sized arrays
ind = sub2ind([K, Nedges], element1, Edge1);

ElmtToElmt = zeros(K,Nedges);
ElmtToEdge = zeros(K,Nedges);
ElmtToElmt(ind) = element2;
ElmtToEdge(ind) = Edge2;
