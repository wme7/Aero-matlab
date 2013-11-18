function [u,Du,eqn,info] = PoissonWG(node,elem,pde,bdFlag,option)
%% POISSONWG Poisson equation: lowest order weak Galerkin element
%
%   u = POISSONWG(node,elem,pde,bdFlag) produces the linear finite element
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
%   u = POISSONWG(node,elem,pde,bdFlag,option) specifies the options.
%    - option.dquadorder: quadrature order for diffusion coefficients
%    - option.fquadorder: quadrature order for computing right hand side f
%    - option.solver
%      'direct': the built in direct solver \ (mldivide)
%      'mg':     multigrid-type solvers mg is used.
%      'amg':    algebraic multigrid method is used.
%      'none':   only assemble the matrix equation but not solve. 
%   The default setting is to use the direct solver for small size problems
%   and multigrid solvers for large size problems. For more options on the
%   multigrid solver mg, type help mg.
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [u,A] = POISSONWG(node,elem,pde,bdFlag) returns also the
%   non-modified stiffness matrix A, which is semi-definite. The kernel of
%   A consists of constant vectors. The matrix A can be used to evulate the
%   bilinear form A(u,v) = u'*A*v, especially the enery norm of a finite
%   element function u by sqrt(u'*A*u).
%
%   [u,A,eqn] = POISSONWG(node,elem,pde,bdFlag) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%   The solution u = AD\b. The output eqn can be used to test other solvers.
%
%   [u,A,eqn,info] = POISSONWG(node,elem,pde,bdFlag) returns also the
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
%     squarePoissonWG;
%
%   See also Poisson, squarePoisson, Lshape, crack, mg
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
NT = size(elem,1);

tic;  % record assembling time
%% Diffusion coefficient
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                                 % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)       % d is a function   
    [lambda,weight] = quadpts(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
        K = K + weight(p)*pde.d(pxy);      
   end
end

%% Construct data structure 
[elem2edge,edge] = dofedge(elem);
NE = size(edge,1); 
Ndof = NT + NE;
elem2dof = NT + elem2edge;

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);

% compute ct2 = 1/mean(||x-xc||^2)
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
ct2 = 3./sum((mid1 - center).^2 + (mid2 - center).^2 + (mid3 - center).^2,2);
[Dphi,area] = gradbasis(node,elem);
clear center mid1 mid2 mid3

% Mbb: edge - edge           
for i = 1:3
    for j = i:3
        % local to global index map
        ii = double(elem2dof(:,i));
        jj = double(elem2dof(:,j));
        % local stiffness matrix
        Aij = 4*dot(Dphi(:,:,i),Dphi(:,:,j),2).*area + 4/9*ct2.*area;
        if ~isempty(pde.d)
            Aij = K.*Aij;
        end        
        if (j==i)
            A = A + sparse(ii,jj,Aij,Ndof,Ndof);
        else
            A = A + sparse([ii,jj],[jj,ii],[Aij; Aij],Ndof,Ndof);        
        end        
    end
end

% Mob: interior - edge
Aij = -4/3*ct2.*area;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
Mob = sparse([(1:NT)', (1:NT)', (1:NT)'], ...
             double(elem2dof(:)), [Aij, Aij, Aij], Ndof, Ndof);
A = A + Mob + Mob';

% Moo: diagonal of interor
Aij = 4*ct2.*area;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
A =  A + sparse(1:NT, 1:NT, Aij, Ndof, Ndof);

clear K Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts(option.fquadorder);
	nQuad = size(lambda,1);
    bt = zeros(NT,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        bt = bt + weight(p)*fp;
    end
    bt = bt.*area;
    b(1:NT) = bt;
end
clear pxy bt

%% Set up boundary conditions
if nargin<=3, bdFlag = []; end
[AD,b,u,freeDof,isPureNeumann] = getbdWG(A,b);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeDof), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if NE <= 2e3  % Direct solver for small size systems
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
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        % eleminate elementwise dof
        option.solver = 'CG';
        option.x0 = u(NT+1:end);
        Aoinv = spdiags(1./diag(AD(1:NT,1:NT)),0,NT,NT);
        Aob = AD(1:NT,NT+1:end);
        Abo = Aob';
        Abb = AD(NT+1:end,NT+1:end);
        Abbm = Abb - Abo*Aoinv*Aob;
        bm = -Abo*Aoinv*b(1:NT) + b(NT+1:end);        
        [ub,info] = mg(Abbm,bm,elem,option,edge);
        u(1:NT) = Aoinv*(b(1:NT) - Aob*ub);
        u(NT+1:end) = ub;        
%         option.x0 = u;
%         [u,info] = mg(AD,b,elem,option,edge);        
    case 'amg'
        option.solver = 'CG';
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    uc = sum(u(1:NT).*area);
    u = u - uc;   % normalization for pure Neumann problem
end

%% Output information
eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof);
info.assembleTime = assembleTime;

%% Compute Du
dudx =  u(elem2dof(:,1)).*Dphi(:,1,1) + u(elem2dof(:,2)).*Dphi(:,1,2) ...
      + u(elem2dof(:,3)).*Dphi(:,1,3);
dudy =  u(elem2dof(:,1)).*Dphi(:,2,1) + u(elem2dof(:,2)).*Dphi(:,2,2) ...
      + u(elem2dof(:,3)).*Dphi(:,2,3);         
Du = -2*[dudx, dudy];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdWG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof,isPureNeumann]= getbdWG(A,b)
    %% GETBDCR Boundary conditions for Poisson equation: WG element.
    
    u =zeros(Ndof,1);
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    idxR = (bdFlag(:) == 3);      % index of Robin edges in bdFlag
    if any(idxR)    
        isRobin = false(NE,1);
        isRobin(elem2edge(idxR)) = true;
        Robin = edge(isRobin,:);  % Robin edges  
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        mid = (node(Robin(:,1),:) + node(Robin(:,2),:))/2;
        ii = NT + find(isRobin);  % for WG: edge dof is after elem dof
        ss = pde.g_R(mid).*edgeLength; % exact for linear g_R
        A = A + sparse(ii,ii,ss,Ndof,Ndof);
    end
    
    % Find Dirichlet boundary nodes: fixedEdge
    fixedEdge = []; freeEdge = [];
    if ~isempty(bdFlag)              % find boundary edges
        idxD = (bdFlag(:) == 1);     % all Dirichlet edges in bdFlag
        isFixedEdge = false(NE,1);
        isFixedEdge(elem2edge(idxD)) = true;  % index of fixed boundary edges
        fixedEdge = find(isFixedEdge);
        freeEdge = find(~isFixedEdge);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given in the input
        s = accumarray(elem2edge(:), 1, [NE 1]);
        fixedEdge = find(s == 1);
        freeEdge = find(s == 2);
    end
    isPureNeumann = false;    
    if isempty(fixedEdge) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedEdge = 1;
        freeEdge = (2:NE)';    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedEdge,fixedEdge)=I, AD(fixedEdge,freeEdge)=0, AD(freeEdge,fixedEdge)=0.
    if ~isempty(fixedEdge)
        bdidx = zeros(Ndof,1); 
        bdidx(NT + fixedEdge) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    
    %% Part 2: Find boundary edges and modify the right hand side b
    % Find boundary edges: Neumann
    Neumann = [];
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag
        isNeumann = elem2edge(idxN | idxR); % index of Neumann and Robin edges
        % since boundary integral is also needed for Robin edges
        Neumann = edge(isNeumann,:);      % Neumann edges        
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        s = accumarray(elem2edge(:), 1, [NE 1]);
        Neumann = edge(s == 1,:);
    end

    % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),1);
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            ge = ge+ weightgN(pp)*gNp;
        end
        ge = ge.*el;
        b(NT+isNeumann) = b(NT+isNeumann) + ge;
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedEdge) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)  % pde.g_D could be a numerical array 
            u(fixedEdge) = pde.g_D(fixedEdge); 
        else % pde.g_D is a function handle
            mid = (node(edge(fixedEdge,1),:) + node(edge(fixedEdge,2),:))/2;
            u(NT+fixedEdge) = pde.g_D(mid);
        end
        b = b - A*u;
        b(NT+fixedEdge) = u(NT+fixedEdge);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    
    freeDof = [(1:NT)'; NT+freeEdge];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of PoissonWG
