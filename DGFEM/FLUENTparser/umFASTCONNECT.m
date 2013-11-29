% (1) compute connectivity arrays for triangles
function [ElmtToElmt, ElmtToEdge] = umCONNECT(ElmtToNode)
Nedges = 3;

% (2) find nber of elements, maximum node number, total local edges
K = size(ElmtToNode,1);
Nnodes = max(max(ElmtToNode));
TotalEdges = Nedges*K;

% (3) list of local Edge to local Node connections
vnum = [[1,2];[2,3];[1,3]];

% (4alt) build sparse representation
spEdgeToNode = zeros(2*K*Nedges, 3);
sk = 1;
sk1 = 1;
for k=1:K
    for edge=1:Nedges
        spEdgeToNode(sk1, 1) = sk; 
        spEdgeToNode(sk1, 2) = ElmtToNode(k,vnum(edge,1)); 
        spEdgeToNode(sk1, 3) = 1;
        sk1 = sk1+1;
        spEdgeToNode(sk1, 1) = sk; 
        spEdgeToNode(sk1, 2) = ElmtToNode(k,vnum(edge,2)); 
        spEdgeToNode(sk1, 3) = 1;
        sk1 = sk1+1;
        sk = sk+1;
    end
end
EdgeToNode = spconvert(spEdgeToNode);
clear spEdgeToNode

spNodeToEdge = zeros(2*K*Nedges, 3);
sk = 1;
sk1 = 1;
for k=1:K
    for edge=1:Nedges
        spNodeToEdge(sk1, 1) = ElmtToNode(k,vnum(edge,1)); 
        spNodeToEdge(sk1, 2) = sk; 
        spNodeToEdge(sk1, 3) = 1;
        sk1 = sk1+1;
        spNodeToEdge(sk1, 1) = ElmtToNode(k,vnum(edge,2)); 
        spNodeToEdge(sk1, 2) =  sk; 
        spNodeToEdge(sk1, 3) = 1;
        sk1 = sk1+1;
        sk = sk+1;
    end
end

NodeToEdge = spconvert(spNodeToEdge);
clear spNodeToEdge

pack

% (5) build global Edge to global Edge  array
EdgeToEdge = EdgeToNode*NodeToEdge; % transpose(EdgeToNode);

clear EdgeToNode
clear NodeToEdge

% (6) remove self connections
%EdgeToEdge = EdgeToEdge - 2*speye(TotalEdges,TotalEdges);

% (7) find complete Edge to Edge connections
[Edges1, Edges2] = find(EdgeToEdge==2);
clear EdgeToEdge

ids = find(Edges1~=Edges2);
Edges1 = Edges1(ids);
Edges2 = Edges2(ids);

clear ids

% (8) convert Edge global number to element and local Edge numbers
element1 = floor( (Edges1-1)/Nedges )  + 1;
Edge1    =   mod( (Edges1-1), Nedges ) + 1;
clear Edges1

element2 = floor( (Edges2-1)/Nedges )  + 1;
Edge2    =   mod( (Edges2-1), Nedges ) + 1;
clear Edges2

% (9) rearrange into K x Nedges sized arrays
ind = sub2ind([K, Nedges], element1, Edge1);

ElmtToElmt = zeros(K,Nedges);
ElmtToEdge = zeros(K,Nedges);
ElmtToElmt(ind) = element2;
ElmtToEdge(ind) = Edge2;
