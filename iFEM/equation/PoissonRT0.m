function [u,sigma,eqn,info] = PoissonRT0(node,elem,pde,bdFlag,option)
%% POISSONRT0 Poisson equation: lowest order RT element.
%
%  [u,sigma] = PoissonRT0(node,elem,pde,bdFlag) produces an approximation of
%  the Poisson equation 
%
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N
%
%  in the mixed formulation:
%
%  ------------------------------------------------------------------------
%  Find (\sigma , u) in H_{g_N,\Gamma_N}(div,\Omega)\times L^2(\Omega) s.t. 
%
%  (d^-1\sigma,\tau) - (div \tau, u)  = <\tau*n,g_D>_{\Gamma_D} 
%  \forall \tau in H_{0,\Gamma_N}(div,\Omega) 
%    - (div \sigma, v)                =  -(f,v)  
%  \forall v in L^2(\Omega) 
%
%  where 
%  H_{g,\Gamma}(div,\Omega) = {\sigma \in H(div,\Omega); \sigma*n = g 
%  on \Gamma \subset \partial\Omega }.
%  ------------------------------------------------------------------------
%
%  The unknown sigma = d*grad(u) is approximated using the lowest order
%  Raviart-Thomas element and u by piecewise constant element (with basis 1).
%
%  [u,sigma] = PoissonRT0(node,elem,pde,bdFlag,option) specifies options
%   - option.solver
%     'direct': the built in direct solver \ (mldivide)
%     'dmg':     multigrid-type solvers mg is used.
%     'uzawapcg': PCG for the Schur complement equation
%     'none': only assemble the matrix equation but not solve
%
%   The default setting is to use the direct solver for small size problems
%   and transforming based multigrid solvers for large size problems. 
%
%  Example
%
%    exampleRT0
%
% Created by Ming Wang. Reorganized by Long Chen. Change basis for u.  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end

%% Diffusion coefficient
if ~isfield(pde,'d'), pde.d = []; end
if isfield(pde,'d') && ~isempty(pde.d)
   if isnumeric(pde.d)
      K = pde.d;                   % d is an array
   else                            % d is a function
      center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
      K = pde.d(center);  % take inverse sequencil.             
   end
else
    K = [];
end

%% Data structure
elemold = elem;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dofedge(elem);
NT = size(elem,1); NE = size(edge,1);
[Dlambda,area,elemSign] = gradbasis(node,elem);

%% Assemble matrix 
Nsigma = NE; Nu = NT; Ndof = Nsigma + Nu;

% M. Mass matrix for RT0 element
M = getmassmatvec(elem2edge,area,Dlambda,'RT0',K);

% B. negative divergence operator
B = icdmat(double(elem2edge),elemSign*[1 -1 1]);

% C. zero matrix.
C = sparse(Nu,Nu);

A = [M B';B C];

%% Assemble right hand side.
fu = zeros(Nu,1);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isempty(pde.f)
	[lambda,weight] = quadpts(option.fquadorder);
	nQuad = size(lambda,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
		fu = fu - fp*weight(p);
    end
    fu = fu.*area;
end
clear fp area
F((Nsigma+1):(Ndof),1) = fu;

%% Boundary Conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,F,bigu,freeDof,isPureNeumannBC] = getbdRT0(F);
eqn = struct('M',AD(1:NE,1:NE),'B',AD(NE+1:end,1:NE),'C',AD(NE+1:end,NE+1:end),...
             'f',F(1:NE),'g',F(NE+1:end),'freeDof',freeDof);

%% Solve the linear system.
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else         % MGCG  solver for large size systems
        option.solver = 'dmg';
    end
elseif strcmp(option.solver,'mg')
    option.solver = 'dmg';    
end
solver = option.solver;
% solve
switch lower(solver);
    case 'direct'
      bigu(freeDof) = AD(freeDof,freeDof)\F(freeDof);
        sigma = bigu(1:NE);
        u = bigu(NE+1:end); info =[];
    case 'none'
        sigma = []; u = []; info =[];        
    case 'dmg'
        [sigma,u] = dmg(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);    
    case 'uzawapcg'
        [sigma,u] = uzawapcg(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);
end
info.solverTime = 0; info.assembleTime = 0; info.itStep = 0;
info.stopErr = 0; info.flag = 0;
if isPureNeumannBC == true % post process for u.
    ubar = sum(u.*area)/sum(area);
    u = u - ubar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,bigu,freeDof,isPureNeumannBC] = getbdRT0(F)
    %% GETBDRT0 Boundary conditions for Poisson equation: RT0 element.
    %
    %  Created by Ming Wang. Improved the check of edgeSign by Long Chen.

    %%
    bigu = zeros(Ndof,1);
    
    %% Boundary conditions
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end

    %% Set up bdFlag
    if isempty(bdFlag) % no bdFlag information
       if ~isempty(pde.g_N) % case: Neumann
           bdFlag = setboundary(node,elem,'Neumann');
       elseif ~isempty(pde.g_D) % case: Dirichlet
           bdFlag = setboundary(node,elem,'Dirichlet');
       end
    end

    %% Find Dirichlet and Neumann dofs 
    if ~isempty(bdFlag)
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isNeumann(elem2edge(bdFlag(:)==2)) = true;
        % Direction of boundary edges may not be the outwards normal
        % direction of the domain. edgeSign is introduced to record this
        % inconsistency.
        edgeSign = ones(NE,1);
        idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first edge is on boundary
        edgeSign(elem2edge(idx,1)) = -1;
        idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second edge is on boundary
        edgeSign(elem2edge(idx,2)) = -1;
        idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% first edge is on boundary
        edgeSign(elem2edge(idx,3)) = -1;
    end
    Dirichlet = edge(isDirichlet,:);
    Neumann = edge(isNeumann,:); 
    isBdDof = false(Ndof,1); 
    isBdDof(isNeumann) = true;   % for mixed method, Neumann edges are fixed
    freeDof = find(~isBdDof);

    %% Dirichlet boundary condition (Neumann BC in mixed form)
    %   We need only modify the rhs on dof associated with Dirichlet
    %   boundary. Compute the int_e \Phi\cdot n g_D on the boundary using
    %   quadrature rules.
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end
    if ~isempty(pde.g_D) && any(isDirichlet) 
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambda,weight] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Dirichlet(:,1),:)+...
                  lambda(ip,2)*node(Dirichlet(:,2),:);               
            F(isDirichlet) = F(isDirichlet) + weight(ip)*pde.g_D(pxy);
        end
        F(isDirichlet) = F(isDirichlet).*edgeSign(isDirichlet);
        % no edge length since the basis of sigma contains it.
    end

    %% Neumann boundary condition (Dirichlet BC in mixed form)
    if ~isempty(pde.g_N) && any(isNeumann)
        % modify the rhs to include Dirichlet boundary condition 
        mid = 1/2*(node(Neumann(:,1),:)+node(Neumann(:,2),:));
        ve = node(Neumann(:,1),:)-node(Neumann(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        if isnumeric(pde.g_N)
            evalg_N = pde.g_N;
        else
            evalg_N = pde.g_N(mid);
        end
        bigu(isNeumann) = edgeLength.*evalg_N;
        if ~isempty(pde.d)
            bigu(isNeumann) = pde.d(mid).*bigu(isNeumann);
        end
        bigu(isNeumann) = bigu(isNeumann).*edgeSign(isNeumann);
        F = F - A*bigu;
        F(isNeumann) = bigu(isNeumann);
    end
    
    %% Pure Neumann boundary condition
    isPureNeumannBC = false;
    if ~any(isDirichlet) && any(isNeumann)
        freeDof = freeDof(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
        isBdDof(end) = true;
        isPureNeumannBC = true;
%         F(end) = 0;
    end

    %% Modify the matrix
    %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
    %  matrix AD by enforcing  |AD(bdNode,bdNode)=I, 
    %  AD(bdNode,FreeNode)=0, AD(FreeNode,bdNode)=0|.
    if any(isBdDof)
       bdidx = zeros(Ndof,1); 
       bdidx(isBdDof) = 1;
       Tbd = spdiags(bdidx,0,Ndof,Ndof);
       T = spdiags(1-bdidx,0,Ndof,Ndof);
       AD = T*A*T + Tbd;
    else
       AD = A;
    end
    end
end