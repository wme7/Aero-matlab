
%% 2-D test
clc; clear all; close all
for k = 8:10
    N = 2^k;
    A = delsq(numgrid('S',N));
    testachol;
end

%% 3-D test
clear all; close all
load oilpump;
for k = 1:3
    [node,elem] = uniformrefine3(node,elem);
    A = assemblematrix3(node,elem);
%     A = A + M;
    [bdNode,bdFace,isBdNode] = findboundary3(elem);
    A = A(~isBdNode,~isBdNode);
    testachol;
%     [node,elem] = uniformrefine3(node,elem);
end