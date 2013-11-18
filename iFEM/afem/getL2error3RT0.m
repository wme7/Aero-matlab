function err = getL2error3RT0(node,elem,sigma,sigmah,markedElem)
%% GETL2ERROR3RT0  L2 norm of RT0 element in 3D.
% 
%  The input sigma can now just be function boundle. sigmah is the 
%  flux through faces. 
%  err = getL2error3RT0(node,elem,sigma,sigmah,markedElem)
%
% Example
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];  % nodes
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7]; % elements
%     elem = label3(node,elem);
%     [node,elem] = uniformbisect3(node,elem);
%     bdFace = setboundary3(node,elem,'Dirichlet','y~=-1','Neumann','y==-1');
%     maxIt = 3;
%     pde = mixBCdata3;
%     err = zeros(maxIt,1); N = zeros(maxIt,1);
%     for i =1:maxIt
%         [node,elem,~,bdFace] = uniformbisect3(node,elem,[],bdFace);
%         barycenter = 1/4.*(node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:)+node(elem(:,4),:));
%         uexa = pde.exactu(barycenter);
%         [u,sigma] = Poisson3RT0(node,elem,pde,bdFace);
%         err(i) = getL2error3(node,elem,pde.exactu,u);
%         N(i) = length(u) + length(sigma) ;
%     end
%     r1 = showrate(N,err,2);
%     legend('||\sigma-\sigma_h||',['N^{' num2str(r1) '}'],...
%           'LOCATION','Best')
%
%% TODO: add jump coefficent and sigma is vector array.
% 
% Created by Ming Wang at Jan 17, 2011, modified m-lint May 15, 2011.
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Construct Data Structure
[elem2dof,dofSign] = dof3RT0(elem);
NT = size(elem,1);% Ndof = max(elem2dof(:)); %N = size(node,1); 
[Dlambda,area] = gradbasis3(node,elem);
locFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];

%% Compute square of the L2 error element-wise
[lambda,w] = quadpts3(3); % quadrature order is 3
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ... 
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    sigmap = sigma(pxy);
    sigmahp = zeros(NT,3);
    for l = 1:4 % for each basis
        i = locFace(l,1); j = locFace(l,2); k = locFace(l,3);
        % phi_k = lambda_iRot_j - lambda_jRot_i;
        sigmahp = sigmahp + repmat(double(dofSign(:,l)).*sigmah(elem2dof(:,l)),1,3)*2.*...
                   (lambda(p,i)*cross(Dlambda(:,:,j),Dlambda(:,:,k)) + ...
                    lambda(p,j)*cross(Dlambda(:,:,k),Dlambda(:,:,i)) + ...    
                    lambda(p,k)*cross(Dlambda(:,:,i),Dlambda(:,:,j)));
    end
    err = err + w(p)*sum((sigmap - sigmahp).^2,2);
end
err = err.*area;
% modify the error
err(isnan(err)) = 0;
if (nargin == 5) && ~isempty(markedElem)
    err = err(markedElem); % L2 err on some marked region
end
err = sqrt(sum(err));
