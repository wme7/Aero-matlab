function [u,sigma,eqn,face] = Poisson3RT0(node,elem,pde,bdFlag,solver,varargin)
%% POISSON3RT0 Poisson equation: lowest order Raviart-Thomas element in 3-D.
%
% Example
%    example3RT0
%
% Created by Ming Wang at Dec 27, 2010.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
[elem,bdFlag] = sortelem3(elem,bdFlag);
[elem2dof,face] = dof3face(elem);
NT = size(elem,1); NF = size(face,1);
localFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3]; 

%% Compute element-wise basis.
[Dlambda,volume,elemSign] = gradbasis3(node,elem);
signedVolume = elemSign.*volume;

%% Assemble matrix 
NdofSigma = NF;
Ndofu = NT;
% Part: M. Mass matrix for RT0 element
M = sparse(NdofSigma,NdofSigma);
for i = 1:4
    for j = i:4 
% Local basis
		% local to global index map and its sign
		ii = double(elem2dof(:,i));
		jj = double(elem2dof(:,j));
        i1 = localFace(i,1); i2 = localFace(i,2); i3 = localFace(i,3); % [i1,i2,i3] is the face opposite to vertex i.
        j1 = localFace(j,1); j2 = localFace(j,2); j3 = localFace(j,3);
		% computation of mass matrix --- (phi_i, phi_j) 
		Mij = 1/20*volume*4.*( ...
              (1+(i1==j1))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                               mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
             +(1+(i1==j2))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                               mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
             +(1+(i1==j3))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
                               mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
             +(1+(i2==j1))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                               mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
             +(1+(i2==j2))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                               mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
             +(1+(i2==j3))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
                               mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
             +(1+(i3==j1))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                               mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
             +(1+(i3==j2))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                               mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
             +(1+(i3==j3))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
                               mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2));
        if (j==i)
            M = M + sparse(ii,jj,Mij,NdofSigma,NdofSigma);
        else
            M = M + sparse([ii;jj],[jj;ii],[Mij; Mij],NdofSigma,NdofSigma);        
        end        
% Local basis
    end
end
clear Mij
% Part: B. divergence operator
B = sparse(Ndofu,NdofSigma);
elem2faceSign = elemSign*[1,-1,1,-1];
for i = 1:4  
 	Bi = 1./signedVolume.*double(elem2faceSign(:,i));
    B = B + sparse((1:Ndofu),double(elem2dof(:,i)),Bi,Ndofu,NdofSigma);
end
C = sparse(Ndofu,Ndofu);
A = [M B';B C];

%% Assemble right hand side.% if ~isempty(bdEdge)
fu = zeros(Ndofu,1);
if ~isempty(pde.f)
	[lambda,weight] = quadpts3(2);
	nQuad = size(lambda,1);
	for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
			+ lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxy).*elemSign;
		fu = fu - fp*weight(p);
	end
end
F = [zeros(NdofSigma,1);fu];
clear lambda weight nQuad fp area

%% Boundary condition
[AD,F,bigu,freeNode,isPureNeumannBC] = getbd3RT0(A,F);
eqn = struct('M',AD(1:NF,1:NF),'B',AD(NF+1:end,1:NF),'C',AD(NF+1:end,NF+1:end),...
             'f',F(1:NF),'g',F(NF+1:end),'freeDof',freeNode);

%% Solve the linear system.
% set solver type
if (nargin<=4)
    if (NT < 1e3)    % Direct solver for small size systems
        solver = 'direct';
    else            % Multigrid-type solver for large size systems
        solver = 'UzawaPCG';
    end
end
if strcmp(solver,'notsolve');
    sigma=[]; u =[];
elseif strcmp(solver,'direct');
    bigu(freeNode) = AD(freeNode,freeNode)\F(freeNode);
    sigma = bigu(1:NF);
    u = bigu(NF+1:end);
    u = u./signedVolume;
else strcmp(solver,'UzawaPCG');
    [sigma,u] =uzawapcg(AD(1:NF,1:NF),AD(NF+1:end,1:NF),F(1:NF),F(NF+1:end));
    u = u./signedVolume;
end
if isPureNeumannBC==true % post process for u.
    barycenter = 1/4.*(node(elem(:,1),:) + node(elem(:,2),:) + ...
                       node(elem(:,3),:)+node(elem(:,4),:));
    uexa = pde.exactu(barycenter);
    u = u + uexa(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [A,F,bigu,freeNode,isPureNeumannBC] = getbd3RT0(A,F)
        %% GETBD3RT0 Boundary conditions for Poisson equation: RT0 in 3D.
        %
        % Created by Ming Wang.
        
        %%
        NdofSigma = size(face,1);
        Ndofu = length(F) - NdofSigma;
        bigu = zeros(NdofSigma+Ndofu,1);
        
        %% Boundary conditions
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'d'), pde.d = []; end
        
        %% Set up bdFlag
        if (isempty(bdFlag)) % no bdFlag information
            if ~isempty(pde.g_N)
                bdFlag = setboundary3(node,elem,'Neumann');
            elseif ~isempty(pde.g_D)
                bdFlag = setboundary3(node,elem,'Dirichlet');
            end
            % case: bdFlag = [], pde.g_D = pde.g_N =[];
            % It is equivalent to homogenous Dirichlet boundary condition
        end
        faceSign = ones(NF,1);
        isDirichlet = false(NF,1);
        isNeumann = false(NF,1);
        if ~isempty(bdFlag)
            % Find Dirichlet and Neumann boundary edges
            isDirichlet(elem2dof(bdFlag(:) == 1)) = true;
            isNeumann(elem2dof(bdFlag(:) == 2)) = true;
            isbdFace(isDirichlet) = true;
            isbdFace(isNeumann) = true;
            bdFaceDof = find(isbdFace);
            bdFaceOutDirec(elem2dof(bdFlag(:,1) ~= 0,1),:) = -Dlambda(bdFlag(:,1) ~= 0,:,1);
            bdFaceOutDirec(elem2dof(bdFlag(:,2) ~= 0,2),:) = -Dlambda(bdFlag(:,2) ~= 0,:,2);
            bdFaceOutDirec(elem2dof(bdFlag(:,3) ~= 0,3),:) = -Dlambda(bdFlag(:,3) ~= 0,:,3);
            bdFaceOutDirec(elem2dof(bdFlag(:,4) ~= 0,4),:) = -Dlambda(bdFlag(:,4) ~= 0,:,4);
            ve1 = node(face(isbdFace,2),:) - node(face(isbdFace,1),:);
            ve2 = node(face(isbdFace,3),:) - node(face(isbdFace,1),:);
            nvface = mycross(ve1,ve2);
            faceSign(bdFaceDof(dot(nvface,bdFaceOutDirec(isbdFace,:),2)<0)) = -1;
        end
        Dirichlet = face(isDirichlet,:);
        Neumann = face(isNeumann,:);
        isbigBdNode = false(NdofSigma+Ndofu,1);
        isbigBdNode(1:NdofSigma) = isNeumann;
        isfreeNode = ~isbigBdNode;
        freeNode = find(isfreeNode);
        
        %% Dirichlet boundary condition (Neumann BC in mixed form)
        %   We need only modify the rhs on dof associated with Dirichlet boundary
        %   Compute the integration of g_D on the boundary Face with middle point
        %   quadrature rule. \Phi\cdot n = |F|.
        %   int_F \Phi\cdot n g_D = g_D(F_barycenter)
        if ~isempty(pde.g_D) && (~isempty(Dirichlet))
            barycenter = 1/3*(node(Dirichlet(:,1),:)+node(Dirichlet(:,2),:)+ ...
                              node(Dirichlet(:,3),:));
            F(isDirichlet) = pde.g_D(barycenter).*faceSign(isDirichlet);
        end
        clear barycenter
        
        %% Neumann boundary condition (Dirichlet BC in mixed form)
        if ~isempty(pde.g_N) && any(isNeumann)
            % modify the rhs to include Dirichlet boundary condition
            barycenter = 1/3*(node(Neumann(:,1),:)+node(Neumann(:,2),:)+node(Neumann(:,3),:));
            ve2 = node(Neumann(:,1),:) - node(Neumann(:,3),:);
            ve3 = node(Neumann(:,2),:) - node(Neumann(:,1),:);
            ve2Crossve3 = mycross(ve2,ve3);
            faceArea = 0.5*sqrt(sum(ve2Crossve3.^2,2));
            bigu(isNeumann) = faceArea.*pde.g_N(barycenter).*faceSign(isNeumann); % 2<1,g_N>
            F = F - A*bigu;
            F(isNeumann) = bigu(isNeumann);
        end
        %% Pure Neumann boundary condition
        isPureNeumannBC = false;
        if ~any(isDirichlet) && any(isNeumann)
            freeNode = freeNode(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
            isbigBdNode(end) = true;
            F(end) = 0;
            isPureNeumannBC = true;
        end
        %% Modify the matrix
        %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
        %  matrix AD by enforcing  |AD(bdNode,bdNode)=I,
        %  AD(bdNode,freeNode)=0, AD(freeNode,bdNode)=0|.
        if ~isempty(isNeumann)
            bdidx = zeros(NdofSigma+Ndofu,1);
            bdidx(isbigBdNode) = 1;
            Tbd = spdiags(bdidx,0,NdofSigma+Ndofu,NdofSigma+Ndofu);
            T = spdiags(1-bdidx,0,NdofSigma+Ndofu,NdofSigma+Ndofu);
            A = T*A*T + Tbd;
        end
    end
end
%% TODO: Write m-lint