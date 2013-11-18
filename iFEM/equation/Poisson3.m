function [u,Du,eqn,info] = Poisson3(node,elem,pde,bdFlag,option,varargin)
%% POISSON3 Poisson equation: P1 linear element in 3-D.
%
%   u = Poisson3(node,elem,pde,bdFlag) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
% 
%   The mesh is given by node and elem and the boundary edge is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%   
%   The function Poisson assembes the matrix equation AD*u = b and solves
%   it by the direct solver (small size <= 2e3) or the multigrid solver
%   (large size > 2e3). The Dirichlet boundary condition is built into the
%   matrix AD and the Neumann boundary condition is build into b.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as the elem array. 
% 
%   u = Poisson3(node,elem,pde,bdFlag,HB,option) specifies the options.
%     - option.solver == 'direct': the built in direct solver \ (mldivide)
%     - option.solver == 'mg':     multigrid-type solvers mg is used.
%     - option.solver == 'notsolve': the solution u = u_D. 
%   The default setting is to use the direct solver for small size problems
%   and multigrid solvers for large size problems. For more options on the
%   multigrid solver mg, type help mg.
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [u,Du] = Poisson3(node,elem,pde,bdFlag,option,HB) returns also the
%   non-modified stiffness matrix A, which is semi-definite. The kernel of
%   A consists of constant vectors. The matrix A can be used to evulate the
%   bilinear form A(u,v) = u'*A*v, especially the enery norm of a finite
%   element function u by sqrt(u'*A*u).
%
%   [u,Du,eqn] = Poisson3(node,elem,pde,bdFlag,HB) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%   The solution u = AD\b. The output eqn can be used to test other solvers.
%
%   [u,Du,eqn,info] = Poisson3(node,elem,pde,bdFlag,HB) returns also the
%   information on the assembeling and solver, which includes:
%     - info.assembleTime: time to assemble the matrix equation
%     - info.solverTime:   time to solve the matrix equation
%     - info.itStep:       number of iteration steps for the mg solver
%     - info.error:        l2 norm of the residual b - A*u
%     - info.flag:         flag for the mg solver.
%       flag = 0: converge within max iteration 
%       flag = 1: iterated maxIt times but did not converge
%       flag = 2: direct solver
%       flag = 3: no solve
%
%   Example
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%       [node,elem] = uniformbisect(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = Poisson(node,elem,pde);
%     figure(1); 
%     showresult(node,elem,u);
%     % Non-homogenous Dirichlet boundary condition
%     pde.f = inline('-4*ones(size(p,1),1)','p');
%     pde.g_D = inline('p(:,1).^2 + p(:,2).^2','p');
%     u = Poisson(node,elem,pde);
%     figure(2); 
%     showresult(node,elem,u);
%     % Homogenous Neumann boundary condition
%     clear pde
%     pde.f = inline('pi^2*cos(pi*p(:,1)).*cos(pi*p(:,2))','p');
%     u = Poisson(node,elem,pde);
%     figure(3);
%     showresult(node,elem,u);
%
%   Example
%     cubePoisson;
%
%   See also Poisson3, squarePoisson, Lshape, crack, mg
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
N = size(node,1); NT = size(elem,1);
Ndof = N;

tic;  % record assembling time
%% Diffusion coefficient
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                           % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d) % d is a function   
    [lambda,weight] = quadpts3(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
        K = K + weight(p)*pde.d(pxyz);           
   end
end

%% Compute geometric quantities and gradient of local basis
[Dphi,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:4
    for j = i:4
        Aij = dot(Dphi(:,:,i),Dphi(:,:,j),2).*volume;
        if ~isempty(pde.d)
            Aij = K.*Aij; 
        end 
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],Ndof,Ndof);        
        end        
    end
end
clear K Aij

%% Assemble right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts3(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,4);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
        for j = 1:4
            bt(:,j) = bt(:,j) + weight(p)*phi(p,j)*fp;
        end
    end
    bt = bt.*repmat(volume,1,4);
    b = accumarray(elem(:),bt(:),[Ndof 1]);
end
clear pxyz

%% Set up boundary conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,b,u,freeNode,isPureNeumann] = getbd3(b);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        if nargin>=6
            HB = varargin{1};
        else
            HB = [];
        end
        [u,info] = mg(AD,b,elem,option,HB); 
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchVolume = accumarray(elem(:),repmat(volume,4,1), [N 1]); 
    uc = sum(u.*patchVolume)/4;
    u = u - uc;   % int u = 0
end

%% Output information
eqn = struct('A',AD,'b',b,'freeNode',freeNode);
info.assembleTime = assembleTime;

%% Compute Du
dudx = u(elem(:,1)).*Dphi(:,1,1)+u(elem(:,2)).*Dphi(:,1,2) ...
      +u(elem(:,3)).*Dphi(:,1,3)+u(elem(:,4)).*Dphi(:,1,4);
dudy = u(elem(:,1)).*Dphi(:,2,1)+u(elem(:,2)).*Dphi(:,2,2) ...
      +u(elem(:,3)).*Dphi(:,2,3)+u(elem(:,4)).*Dphi(:,2,4);
dudz = u(elem(:,1)).*Dphi(:,3,1)+u(elem(:,2)).*Dphi(:,3,2) ...
      +u(elem(:,3)).*Dphi(:,3,3)+u(elem(:,4)).*Dphi(:,3,4);
Du = [dudx, dudy, dudz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode,isPureNeumann] = getbd3(b)
    %% Set up of boundary conditions.
    %
    % 1) Modify the matrix for Dirichlet boundary nodes, which are not degree
    % of freedom. Values at these nodes are evaluatation of pde.g_D. The
    % original stiffness matrix A is turn into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    %
    % 2) Modify the right hand side b. The Neumann boundary integral is added
    % to b. For Dirichlet boundary ndoes, b(fixedNode) is the evaluation of
    % pde.g_D.
    %
    % Special attentation should be given for the pure Neumann boundary
    % condition. To enforce the compatible condition, the vector b should have
    % mean value zero. To avoid a singular matrix, the 1st node is chosen as
    % fixedNode. 
    %
    % The order of assigning Neumann and Dirichlet boundary condition is
    % important to get the right setting at the intersection nodes of Dirichlet
    % and Neumann boundary edges.
    %
    % Reference: Long Chen. Finite Element Methods and its Programming. Lecture
    % Notes.

    u = zeros(Ndof,1); 
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    isRobin = (bdFlag(:) == 3);
    if any(isRobin)
        allFace = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
        Robin = allFace(isRobin,:);
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        v12 = node(Robin(:,2),:)-node(Robin(:,1),:);
        v13 = node(Robin(:,3),:)-node(Robin(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        if ~isfield(option,'gRquadorder')
            option.gRquadorder = 2;   % default order exact for linear gR
        end
        [lambdagR,weightgR] = quadpts(option.gRquadorder);
        nQuadgR = size(lambdagR,1);
        bdphi = lambdagR;
        NR = size(Robin,1);
        ss = zeros(NR,3,3);
        % int g_R phi_i phi_j
        for pp = 1:nQuadgR
            ppxyz = lambdagR(pp,1)*node(Robin(:,1),:) ...
                  + lambdagR(pp,2)*node(Robin(:,2),:) ...
                  + lambdagR(pp,3)*node(Robin(:,3),:);
            gRp = pde.g_R(ppxyz);
            for iR = 1:3
                for jR = iR:3
                    ss(:,iR,jR) = ss(:,iR,jR) + ...
                    weightgR(pp)*gRp*bdphi(pp,iR).*bdphi(pp,jR);
                end
            end
        end
        ss(:) = ss(:).*repmat(area,9,1);       
        index = 0;
        for iR = 1:3
            for jR = 1:3
                iiR(index+1:index+NR) = double(Robin(:,iR)); 
                jjR(index+1:index+NR) = double(Robin(:,jR)); 
                if jR>=iR
                    ssR(index+1:index+NR) = ss(:,iR,jR);
                else
                    ssR(index+1:index+NR) = ss(:,jR,iR);
                end
                index = index + NR;
            end
        end
        A = A + sparse(iiR,jjR,ssR,Ndof,Ndof);        
    end

    % Find Dirichlet boundary nodes: fixedNode
    fixedNode = []; freeNode = [];
    if ~isempty(bdFlag) % bdFlag specifies different bd conditions
        [fixedNode,bdFace,isBdNode] = findboundary3(elem,bdFlag);
        freeNode = find(~isBdNode);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given
        [fixedNode,bdFace,isBdNode] = findboundary3(elem);
        freeNode = find(~isBdNode);
    end
    isPureNeumann = false;
    if isempty(fixedNode) && isempty(Robin) % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedNode = 1;
        freeNode = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    if ~isempty(fixedNode)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end        

    %% Part 2: Find boundary faces and modify the load b
    % Find boundary faces: Neumann
    Neumann = [];
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        Neumann = bdFace;       
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        [tempvar,Neumann] = findboundary(elem); %#ok<ASGLU> %TODO: is this findboundary3?
    end

    % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),3);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                  + lambdagN(pp,2)*node(Neumann(:,2),:) ...
                  + lambdagN(pp,3)*node(Neumann(:,3),:);
            gNp = pde.g_N(ppxyz);
            for iN = 1:3
                ge(:,iN) = ge(:,iN) + weightgN(pp)*phigN(pp,iN)*gNp;
            end
        end
        ge = ge.*repmat(area,1,3);
        b = b + accumarray(Neumann(:),ge(:),[Ndof,1]); 
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedNode) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)
            u(fixedNode) = pde.g_D(fixedNode);
        else
            u(fixedNode) = pde.g_D(node(fixedNode,:));
        end
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedNode) = u(fixedNode);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;          % 1 is fixedNode and set u(1) = 0
    end
    end % end of getbd3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Poisson3