function [ElmtToElmt, ElmtToEdge] = umSLOWCONNECT(ElmtToNode)

% input: ElmtToNode  assumed to be K x 3, i.e. K rows of 3 vertex numbers each
K = size(ElmtToNode,1); % find number of triangles

% output matrices:
ElmtToElmt = zeros(K, 3);
ElmtToEdge = zeros(K,3);

% local vertices on each local edge
vnum = [[1,2];[2,3];[3,1]];

% nested loop approach to find element to element connectivity
for elmt1=1:K
    for edge1=1:3
        % find the two vertices on this
        v1a = ElmtToNode(elmt1,vnum(edge1,1));
        v1b = ElmtToNode(elmt1,vnum(edge1,2));
        for elmt2=1:K
            if(elmt1 ~= elmt2) % make sure we are looking at a neighbor
                for edge2=1:3
                    v2a = ElmtToNode(elmt2,vnum(edge2,1));
                    v2b = ElmtToNode(elmt2,vnum(edge2,2));
            
                    % just need to check one case
                    % (since all elements are counter clockwise)
                    if( (v1a == v2b) & (v1b == v2a))
                        ElmtToElmt(elmt1, edge1) = elmt2;
                        ElmtToEdge(elmt1, edge1) = edge2;
                    end
                end
            end
        end
    end
end
