function [u,sigma,eqn,info] = PoissonBDM1(node,elem,pde,bdFlag,option)
%% POISSONBDM1 Poisson equation: linear BDM element.
%
%  [u,sigma] = PoissonBDM1(node,elem,pde,bdFlag) produces an approximation of
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
%  (d^-1\sigma,\tau) + (div \tau, u)  = <\tau*n,g_D>_{\Gamma_D} 
%  \forall \tau in H_{0,\Gamma_N}(div,\Omega) 
%  (div \sigma, v)                =  -(f,v)  
%  \forall v in L^2(\Omega) 
%
%  where 
%  H_{g,\Gamma}(div,\Omega) = {\sigma \in H(div,\Omega); \sigma*n = g 
%  on \Gamma \subset \partial\Omega }
%  ------------------------------------------------------------------------
%
%  The unknown sigma = d*grad(u) is approximated using the Brezzi-Douglas-Michel 
%  element and u  by piecewise constant element (with basis 1).
%
%  [u,sigma] = PoissonBDM1(node,elem,pde,bdFlag,option) specifies options
%   - option.solver
%     'direct': the built in direct solver \ (mldivide)
%     'dmg':     multigrid-type solvers mg is used.
%     'uzawapcg': PCG for the Schur complement equation
%     'none': only assemble the matrix equation but not solve
%
%   The default setting is to use the direct solver for small size problems
%   and transforming based multigrid solvers for large size problems. 
%
% Example
%
%    exampleBDM1
%
% Created by Ming Wang. Revised by Long Chen.
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
[elem,bdFlag] = sortelem(elem,bdFlag);  % use ascend ordering
[elem2edge,edge] = dofedge(elem);
NT = size(elem,1); NE = size(edge,1);

%% Local basis
[Dlambda,area,elemSign] = gradbasis(node,elem);

%% Assemble matrix 
Nsigma = 2*NE; Nu = NT; Ndof = Nsigma + Nu;

% Part: M. Mass matrix for BDM1 element
M = getmassmatvec(elem2edge,area,Dlambda,'BDM1',K);

% Part: B. divergence operator 
B = icdmat(double(elem2edge),elemSign*[1 -1 1]);
B = [B sparse(NT,NE)];

% Part: C. zero matrix.
C = sparse(Nu,Nu);

A = [M B';B C];

%% Assemble right hand side.
fu = zeros(Nu,1);
if isfield(pde,'f') && isreal(pde.f) && (pde.f==0)
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
F(Nsigma+1:(Nsigma+Nu),1) = fu;

%% Boundary conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,F,bigu,freeDof,isPureNeumannBC] = getbdBDM1(F);
eqn = struct('M',AD(1:2*NE,1:2*NE),'B',AD(2*NE+1:end,1:2*NE),'C',AD(2*NE+1:end,2*NE+1:end),...
             'f',F(1:2*NE),'g',F(2*NE+1:end),'freeDof',freeDof);

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
        sigma = bigu(1:Nsigma);
        u = bigu(Nsigma+1:end); info =[];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,bigu,freeDof,isPureNeumannBC] = getbdBDM1(F)
        %% GETBDBDM1 Boundary conditions for Poisson equation: BDM1 element.
        %
        %  Created by Ming Wang.
        
        %%
        bigu = zeros(Nsigma+Nu,1);
        
        %% Boundary conditions
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        
        %% Set up bdFlag
        if (isempty(bdFlag)) % no bdFlag information
            if ~isempty(pde.g_N) % case:
                bdFlag = setboundary(node,elem,'Neumann');
            elseif ~isempty(pde.g_D) % case:
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
        idpsi = NE+(1:NE);
        isBdDof = false(Nsigma+Nu,1);
        isBdDof(isNeumann) = true;
        isBdDof(idpsi(isNeumann)) = true;
        freeDof = find(~isBdDof);
       
        %% Dirichlet boundary condition (Neumann BC in mixed form)
        %   We need only modify the rhs on dof associated with Dirichlet boundary
        %   Compute the integration of g_D on the boundary with simpson
        %   quadrature rule.
        if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
            pde.g_D = [];
        end
        if ~isempty(pde.g_D) && (any(isDirichlet))
            [lambda,weight] = quadpts1(4);
            nQuad = size(lambda,1);
            % <\phi \cdot n, g_D> = 1/|e_{i,j}|\int e_{i,j} g_D ds
            % <\psi\cdot n,g_D>=1/|e_{i,j}|\int e_{i,j}(\lambda_j-\lambda_i)g_Dds
            for ip = 1:nQuad
                pxy = lambda(ip,1)*node(Dirichlet(:,1),:)+lambda(ip,2)*node(Dirichlet(:,2),:);
                F(isDirichlet) = F(isDirichlet) + weight(ip)*pde.g_D(pxy);
                F(idpsi(isDirichlet)) = F(idpsi(isDirichlet)) + ...
                                        weight(ip)*(lambda(ip,1)-lambda(ip,2))*pde.g_D(pxy);
            end
            F(isDirichlet) = F(isDirichlet).*edgeSign(isDirichlet);
            F(idpsi(isDirichlet)) = F(idpsi(isDirichlet)).*edgeSign(isDirichlet);
        end
        %% Neumann boundary condition (Dirichlet BC in mixed form)
        if ~isempty(pde.g_N) && any(isNeumann)
            % modify the rhs to include Dirichlet boundary condition
            % Here we use simpson formula.
            % The dual functional for bases phi and psi are:
            % \Phi(u) = \int_e{i,j} u \cdot n_{i,j} ds
            % \Psi(u) = 3\int_e{i,j} u \cdot n_{i,j}(\lambda_i-\lambda_j) ds
            ve = node(Neumann(:,1),:)-node(Neumann(:,2),:);
            edgeLength = sqrt(sum(ve.^2,2));
            [lambda,weight] = quadpts1(4);
            nQuad = size(lambda,1);
            bigu(isNeumann) = 0;bigu(idpsi(isNeumann))=0;
            for ip = 1:nQuad
                pxy = lambda(ip,1)*node(Neumann(:,1),:) + lambda(ip,2)*node(Neumann(:,2),:);
                bigu(isNeumann) = bigu(isNeumann) + weight(ip)*pde.g_N(pxy).*edgeLength;
                bigu(idpsi(isNeumann)) = bigu(idpsi(isNeumann)) + ...
                                         weight(ip)*3*(lambda(ip,1)-lambda(ip,2))*pde.g_N(pxy).*edgeLength;
            end
            bigu(isNeumann) = bigu(isNeumann).*edgeSign(isNeumann);
            bigu(idpsi(isNeumann)) = bigu(idpsi(isNeumann)).*edgeSign(isNeumann);
            F = F - A*bigu;
            F(isNeumann) = bigu(isNeumann); F(idpsi(isNeumann))=bigu(idpsi(isNeumann));
        end
        %% Pure Neumann boundary condition
        isPureNeumannBC = false;
        if ~any(isDirichlet) && any(isNeumann)
            freeDof = freeDof(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
            isBdDof(end) = true;
            isPureNeumannBC = true;
        end
        %% Modify the matrix
        %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
        %  matrix AD by enforcing  |AD(bdNode,bdNode)=I,
        %  AD(bdNode,freeDof)=0, AD(freeDof,bdNode)=0|.
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
