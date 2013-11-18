function T = auxstructurec(elem)
%% AUXSTRUCTUREC auxiliary structure for a 2-D triangulation (C-code).
%
%  auxstructurec does the same as auxstructure except the core part is
%  written in C. To use it, please complie auxstructureccode by type
%  
%  	mex auxstructureccode
%
%  See also auxstructure, auxstructure3.
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

[neighbor,elem2edge,edge,bdEdge] = auxstructureccode(int32(elem));
edge = sort(edge,2);
NT = size(elem,1);
NE = size(edge,1);
edge2elem = zeros(NE,4);
allEdgeIdx = elem2edge(:);
allLEdgeIdx = [ones(NT,1);2*ones(NT,1);3*ones(NT,1)];
edge2elem(allEdgeIdx,1) = [1:NT,1:NT,1:NT]';
edge2elem(allEdgeIdx,3) = allLEdgeIdx;
edge2elem(allEdgeIdx(end:-1:1),2) = [NT:-1:1,NT:-1:1,NT:-1:1]';
edge2elem(allEdgeIdx(end:-1:1),4) = allLEdgeIdx(end:-1:1);
T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,...
           'edge2elem', edge2elem,'bdEdge',bdEdge);