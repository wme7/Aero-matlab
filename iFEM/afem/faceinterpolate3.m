function fI = faceinterpolate3(f,node,face,quadOrder)
%% FACEINTERPOLATE3 interpolate to face elements RT0.
%
% uI = faceinterpolate3(u,node,face) interpolates a given function u
% into the lowesr order RT0. The coefficient is given by the face integral 
% int_f u*n ds. 
%
% Example
%
%   
%    [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%    maxIt = 3;
%    pde = polynomialdata0;
%    err = zeros(maxIt,1); 
%    N = zeros(maxIt,1);
%    for i =1:maxIt
%      [node,elem] = uniformrefine3(node,elem);
%      uI = faceinterpolate3(pde.u,node,elem);
%      err(i) = getL2error3RT0(node,elem,pde.u,uI);
%      N(i) = size(u,1);
%    end
%     figure;
%     showrate(N,err,2,'r-+','||u - u_I||');
%
% See also edgeinterpolate, edgeinterpolate1, edgeinterpolate2
%
% Created by Lin Zhong and Ming Wang at July, 2012.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%%  
NF = size(face,1);
v12 = node(face(:,2),:) - node(face(:,1),:);
v13 = node(face(:,3),:) - node(face(:,1),:);
nVec = mycross(v12,v13);

%% assemble the right hand side
if ~exist('quadOrder','var'), quadOrder = 2; end
[lambda,weight] = quadpts(quadOrder);
nQuad = size(lambda,1);
bt = zeros(NF,3);
for p = 1:nQuad
    pxy = lambda(p,1)*node(face(:,1),:) ...
		+ lambda(p,2)*node(face(:,2),:) ...
		+ lambda(p,3)*node(face(:,3),:);
    bt = bt + weight(p)*f(pxy);
end
fI = dot(nVec,bt,2)/2;