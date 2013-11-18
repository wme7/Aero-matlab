function [isM,As] = mis(A,theta)
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

if ~exist('theta','var'), theta = 0.025; end

N = size(A,1);
%% Generate strong connectness matrix
D = spdiags(1./sqrt(diag(A)),0,N,N);
Am = D*A*D;
[im,jm,sm] = find(Am);
idx = (-sm > theta);    % delete weakly connect off-diagonal and diagonal
As = sparse(im(idx),jm(idx),sm(idx),N,N); % matrix for strong connectness
%%
% The diagonal of Am is 1. The negative off-diagonal measures the
% diffusivity. The positive off-diagonal is filtered out. 

%% Compute degree of vertex
deg = sum(spones(As));  % number of strongly connected neighbors
deg = full(deg');
idx = (deg>0);
deg(idx) = deg(idx) + 0.1*rand(sum(idx),1); % break the equal degree case

%% Find an approximate maximal independent set and put to F set
isM = false(N,1);       % M: marked nodes
isM(~idx) = true;       % mark all isolated nodes
% isLeft = false(N,1);    % F: fine node
isU = ~isM;        % U: undecided set
while sum(isU) > sqrt(N)
    % Mark all undecided nodes
    isS = isU;    % S: selected set, changing in the procedure
    S = find(isU); 

    % Find marked nodes with local minimum degree
    [i,j] = find(triu(As(S,S),1));  % i,j and i<j: edges of subgraph S
    idx = deg(S(i)) <= deg(S(j));   % compare degree of connected vertices
    isS(S(j(idx))) = false;         % remove vertices with smaller degree 
    isS(S(i(~idx))) = false;        % if degrees are equal, keep the nodes with smaller index
    isM(isS) = true;                % mark selected nodes

    % Remove selected nodes and neighboring nodes from undecided set
    [i,j] = find(As(:,isM)); %#ok<*NASGU>
    isU(i) = false;     
    isU(isM) = false;
end
% fprintf('Number of marked Nodes: %6.0u\n',sum(isM));